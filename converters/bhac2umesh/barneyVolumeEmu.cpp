// barneyVolumeEmu: CPU emulation of barney's unstructured-mesh volume path
// (diagnostics only). For a set of .umesh parts (one per data rank) and a
// haystack .xf transfer function it reproduces, in float arithmetic like the GPU:
//   - UMeshField::buildInitialMacroCells (MC grid dims/origin/spacing, element
//     bounding boxes rasterised with rasterBox)          native/umesh/common/UMeshField.cu
//   - mapMCs / TransferFunction::DD::majorant           native/volume/MCGrid.cu, TransferFunction.h
//   - MCVolumeAccel::isProg: boxTest, dda3 over the MC grid, Woodcock::sampleRange
//     with fastLog                                        native/volume/MCAccelerator.h, DDA.h, Volume.h
//   - hex/prism point location + interpolation (Newton, VTK order)
//                                                        native/umesh/common/ElementIntersection.h
// and renders the first volume "hit" of primary rays (closest over all parts,
// like the data-parallel merge). No lighting: the pixel is the TF colour at the
// hit (averaged over spp), plus a hit-depth image. Mode --global-majorant uses
// one majorant per part instead of the MC grid (an unbiased reference for the
// same sampler), --exact-log uses logf instead of fastLog.
//
// usage: barneyVolumeEmu -xf tf.xf -o out_prefix [--camera px py pz lx ly lz ux uy uz]
//        [-fovy 40] [-res W H] [-spp N] [--global-majorant] [--exact-log] part0.umesh ...

#include "umesh/UMesh.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using umesh::UMesh;

namespace {

  struct V3 { float x, y, z; };
  inline V3 operator+(V3 a, V3 b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
  inline V3 operator-(V3 a, V3 b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
  inline V3 operator*(V3 a, float s) { return { a.x * s, a.y * s, a.z * s }; }
  inline float dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
  inline V3 cross(V3 a, V3 b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }
  inline float det(V3 a, V3 b, V3 c) { return dot(a, cross(b, c)); }
  inline float get(V3 v, int d) { return d == 0 ? v.x : (d == 1 ? v.y : v.z); }
  inline void set(V3 &v, int d, float f) { if (d == 0) v.x = f; else if (d == 1) v.y = f; else v.z = f; }

  // ---------------- transfer function (haystack .xf + barney TF) ----------------
  struct XF {
    std::vector<float> r, g, b, a;
    float lo = 0, hi = 1, baseDensity = 1;
    void load(const std::string &fn)
    {
      std::ifstream in(fn, std::ios::binary);
      if (!in) throw std::runtime_error("cannot open " + fn);
      size_t magic; float bd; float absD[2], relD[2]; int n;
      in.read((char *)&magic, 8); in.read((char *)&bd, 4);
      in.read((char *)absD, 8); in.read((char *)relD, 8); in.read((char *)&n, 4);
      std::vector<float> cm(4 * (size_t)n);
      in.read((char *)cm.data(), cm.size() * 4);
      for (int i = 0; i < n; i++) { r.push_back(cm[4*i]); g.push_back(cm[4*i+1]); b.push_back(cm[4*i+2]); a.push_back(cm[4*i+3]); }
      const float s = absD[1] - absD[0];
      lo = absD[0] + relD[0] / 100.f * s; hi = absD[0] + relD[1] / 100.f * s;
      // haystack: unitDistance = 1.05^(bd-100); barney anari: densityScale = 1/unitDistance
      baseDensity = 1.f / powf(1.05f, bd - 100.f);
      printf("xf %s: %d entries, domain [%g, %g], baseDensity(xf) %g -> density scale %g\n",
             fn.c_str(), n, lo, hi, bd, baseDensity);
    }
    int n() const { return (int)a.size(); }
    // TransferFunction::DD::map
    void map(float s, float *rgb, float &w) const
    {
      float f = (s - lo) / (hi - lo);
      f = std::min(1.f, std::max(0.f, f)) * (n() - 1);
      int idx = std::min(n() - 2, std::max(0, (int)f));
      f -= idx;
      rgb[0] = (1 - f) * r[idx] + f * r[idx + 1];
      rgb[1] = (1 - f) * g[idx] + f * g[idx + 1];
      rgb[2] = (1 - f) * b[idx] + f * b[idx + 1];
      w = ((1 - f) * a[idx] + f * a[idx + 1]) * baseDensity;
    }
    // TransferFunction::DD::majorant
    float majorant(float rl, float ru) const
    {
      if (rl > ru) return 0.f;
      float fl = std::min(1.f, std::max(0.f, (rl - lo) / (hi - lo))) * (n() - 1);
      float fh = std::min(1.f, std::max(0.f, (ru - lo) / (hi - lo))) * (n() - 1);
      int il = std::min(n() - 2, std::max(0, (int)fl)), ih = std::min(n() - 2, std::max(0, (int)fh));
      fl -= il; fh -= ih;
      float m = std::max((1 - fl) * a[il] + fl * a[il + 1], (1 - fh) * a[ih] + fh * a[ih + 1]);
      for (int i = il + 1; i <= ih; i++) m = std::max(m, a[i]);
      return m * baseDensity;
    }
  };

  // ---------------- element interpolation (ElementIntersection.h) ----------------
  void hexW(const float *p, float *sf, float *d)
  {
    const float rm = 1 - p[0], sm = 1 - p[1], tm = 1 - p[2];
    sf[0] = rm*sm*tm; sf[1] = p[0]*sm*tm; sf[2] = p[0]*p[1]*tm; sf[3] = rm*p[1]*tm;
    sf[4] = rm*sm*p[2]; sf[5] = p[0]*sm*p[2]; sf[6] = p[0]*p[1]*p[2]; sf[7] = rm*p[1]*p[2];
    d[0] = -sm*tm; d[1] = sm*tm; d[2] = p[1]*tm; d[3] = -p[1]*tm;
    d[4] = -sm*p[2]; d[5] = sm*p[2]; d[6] = p[1]*p[2]; d[7] = -p[1]*p[2];
    d[8] = -rm*tm; d[9] = -p[0]*tm; d[10] = p[0]*tm; d[11] = rm*tm;
    d[12] = -rm*p[2]; d[13] = -p[0]*p[2]; d[14] = p[0]*p[2]; d[15] = rm*p[2];
    d[16] = -rm*sm; d[17] = -p[0]*sm; d[18] = -p[0]*p[1]; d[19] = -rm*p[1];
    d[20] = rm*sm; d[21] = p[0]*sm; d[22] = p[0]*p[1]; d[23] = rm*p[1];
  }
  void prismW(const float *p, float *sf, float *d)
  {
    sf[0] = (1 - p[0] - p[1]) * (1 - p[2]); sf[1] = p[0] * (1 - p[2]); sf[2] = p[1] * (1 - p[2]);
    sf[3] = (1 - p[0] - p[1]) * p[2]; sf[4] = p[0] * p[2]; sf[5] = p[1] * p[2];
    d[0] = -1 + p[2]; d[1] = 1 - p[2]; d[2] = 0; d[3] = -p[2]; d[4] = p[2]; d[5] = 0;
    d[6] = -1 + p[2]; d[7] = 0; d[8] = 1 - p[2]; d[9] = -p[2]; d[10] = 0; d[11] = p[2];
    d[12] = -1 + p[0] + p[1]; d[13] = -p[0]; d[14] = -p[1]; d[15] = 1 - p[0] - p[1]; d[16] = p[0]; d[17] = p[1];
  }
  bool eltValue(int n, const V3 *V, const float *S, V3 P, float &val)
  {
    float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
    for (int i = 0; i < n; i++) for (int d = 0; d < 3; d++) { lo[d] = std::min(lo[d], get(V[i], d)); hi[d] = std::max(hi[d], get(V[i], d)); }
    const float diag = std::sqrt((hi[0]-lo[0])*(hi[0]-lo[0]) + (hi[1]-lo[1])*(hi[1]-lo[1]) + (hi[2]-lo[2])*(hi[2]-lo[2]));
    const float detTol = (n == 8) ? 1e-10f * diag : 1e-6f;
    float p[3] = { .5f, .5f, .5f }, sf[8], d[24];
    bool conv = false;
    for (int it = 0; it < 10; it++) {
      if (n == 8) hexW(p, sf, d); else prismW(p, sf, d);
      V3 f = { 0, 0, 0 }, r = { 0, 0, 0 }, s = { 0, 0, 0 }, t = { 0, 0, 0 };
      for (int i = 0; i < n; i++) { f = f + V[i] * sf[i]; r = r + V[i] * d[i]; s = s + V[i] * d[i + n]; t = t + V[i] * d[i + 2 * n]; }
      f = f - P;
      const float dd = det(r, s, t);
      if (std::fabs(dd) < detTol) return false;
      const float d0 = det(f, s, t) / dd, d1 = det(r, f, t) / dd, d2 = det(r, s, f) / dd;
      p[0] -= d0; p[1] -= d1; p[2] -= d2;
      if (std::fabs(d0) < 1e-4f && std::fabs(d1) < 1e-4f && std::fabs(d2) < 1e-4f) { conv = true; break; }
      if (std::fabs(p[0]) > 1e6f || std::fabs(p[1]) > 1e6f || std::fabs(p[2]) > 1e6f) return false;
    }
    if (!conv) return false;
    const float L = -1e-6f, H = 1 + 1e-6f;
    bool in = p[0] >= L && p[0] <= H && p[1] >= L && p[1] <= H && p[2] >= L && p[2] <= H;
    if (n == 6) in = in && p[0] + p[1] <= H;
    if (!in) return false;
    if (n == 8) hexW(p, sf, d); else prismW(p, sf, d);
    val = 0;
    for (int i = 0; i < n; i++) val += sf[i] * S[i];
    return true;
  }

  // ---------------- one data rank ----------------
  struct Part {
    std::vector<V3> vtx;
    std::vector<float> scal;
    std::vector<int> idx, ofs;      // element indices, offsets
    std::vector<uint8_t> nv;        // 8 or 6
    V3 lo, hi;
    // MC grid
    int dims[3];
    V3 origin, spacing;
    std::vector<float> majorant;
    float globalMajorant = 0;
    // point-location grid (CPU stand-in for cuBQL; only finds candidates)
    int gd[3];
    V3 gsp;
    std::vector<uint32_t> gstart, gitems;

    void elementBox(int e, float *blo, float *bhi, float &smin, float &smax) const
    {
      for (int d = 0; d < 3; d++) { blo[d] = 1e30f; bhi[d] = -1e30f; }
      smin = 1e30f; smax = -1e30f;
      for (int i = 0; i < nv[e]; i++) {
        const int v = idx[ofs[e] + i];
        for (int d = 0; d < 3; d++) { blo[d] = std::min(blo[d], get(vtx[v], d)); bhi[d] = std::max(bhi[d], get(vtx[v], d)); }
        smin = std::min(smin, scal[v]); smax = std::max(smax, scal[v]);
      }
    }

    void load(const std::string &fn, const XF &xf)
    {
      UMesh::SP m = UMesh::loadFrom(fn);
      for (auto &v : m->vertices) vtx.push_back({ v.x, v.y, v.z });
      scal = m->perVertex->values;
      // same order as haystack AnariDeviceRenderer::create(UMesh): tets, pyrs, wedges, hexes
      for (auto &w : m->wedges) { ofs.push_back((int)idx.size()); nv.push_back(6); for (int i = 0; i < 6; i++) idx.push_back(w[i]); }
      for (auto &h : m->hexes) { ofs.push_back((int)idx.size()); nv.push_back(8); for (int i = 0; i < 8; i++) idx.push_back(h[i]); }
      if (!m->tets.empty() || !m->pyrs.empty()) throw std::runtime_error("tets/pyrs not emulated");
      const int ne = (int)nv.size();
      lo = { 1e30f, 1e30f, 1e30f }; hi = { -1e30f, -1e30f, -1e30f };
      for (int e = 0; e < ne; e++) {
        float a[3], b[3], s0, s1;
        elementBox(e, a, b, s0, s1);
        for (int d = 0; d < 3; d++) { set(lo, d, std::min(get(lo, d), a[d])); set(hi, d, std::max(get(hi, d), b[d])); }
      }
      // UMeshField::buildInitialMacroCells
      const V3 size = hi - lo;
      const float maxWidth = std::max(size.x, std::max(size.y, size.z));
      const int MC = 400 + int(sqrtf(ne / 1000.f));
      for (int d = 0; d < 3; d++) dims[d] = 1 + int(get(size, d) * ((MC - 1) / maxWidth));
      origin = lo;
      spacing = { size.x / dims[0], size.y / dims[1], size.z / dims[2] };
      const size_t nmc = (size_t)dims[0] * dims[1] * dims[2];
      std::vector<float> rlo(nmc, 1e30f), rhi(nmc, -1e30f);
      for (int e = 0; e < ne; e++) {
        float a[3], b[3], s0, s1;
        elementBox(e, a, b, s0, s1);
        int l[3], h[3];
        for (int d = 0; d < 3; d++) {   // rasterBox
          l[d] = int((a[d] - get(origin, d)) * (1.f / get(spacing, d)));
          h[d] = int((b[d] - get(origin, d)) * (1.f / get(spacing, d)));
          l[d] = std::min(std::max(l[d], 0), dims[d] - 1);
          h[d] = std::min(std::max(h[d], 0), dims[d] - 1);
        }
        for (int z = l[2]; z <= h[2]; z++) for (int y = l[1]; y <= h[1]; y++) for (int x = l[0]; x <= h[0]; x++) {
          const size_t c = x + (size_t)dims[0] * (y + (size_t)dims[1] * z);
          rlo[c] = std::min(rlo[c], s0); rhi[c] = std::max(rhi[c], s1);
        }
      }
      majorant.resize(nmc);
      for (size_t c = 0; c < nmc; c++) { majorant[c] = xf.majorant(rlo[c], rhi[c]); globalMajorant = std::max(globalMajorant, majorant[c]); }
      // candidate grid
      const int G = 48;
      for (int d = 0; d < 3; d++) gd[d] = std::max(1, int(G * get(size, d) / maxWidth));
      gsp = { size.x / gd[0], size.y / gd[1], size.z / gd[2] };
      const size_t ng = (size_t)gd[0] * gd[1] * gd[2];
      std::vector<std::vector<uint32_t>> lists(ng);
      for (int e = 0; e < ne; e++) {
        float a[3], b[3], s0, s1;
        elementBox(e, a, b, s0, s1);
        int l[3], h[3];
        for (int d = 0; d < 3; d++) {
          l[d] = std::min(gd[d] - 1, std::max(0, int((a[d] - get(lo, d)) / get(gsp, d))));
          h[d] = std::min(gd[d] - 1, std::max(0, int((b[d] - get(lo, d)) / get(gsp, d))));
        }
        for (int z = l[2]; z <= h[2]; z++) for (int y = l[1]; y <= h[1]; y++) for (int x = l[0]; x <= h[0]; x++)
          lists[x + (size_t)gd[0] * (y + (size_t)gd[1] * z)].push_back(e);
      }
      gstart.assign(ng + 1, 0);
      for (size_t c = 0; c < ng; c++) gstart[c + 1] = gstart[c] + (uint32_t)lists[c].size();
      gitems.reserve(gstart[ng]);
      for (auto &l : lists) gitems.insert(gitems.end(), l.begin(), l.end());
    }

    // UMeshCuBQLSampler::DD::sample (NaN outside)
    float sample(V3 P) const
    {
      int c[3];
      for (int d = 0; d < 3; d++) {
        const float u = (get(P, d) - get(lo, d)) / get(gsp, d);
        if (u < -1e-3f || u > gd[d] + 1e-3f) return NAN;
        c[d] = std::min(gd[d] - 1, std::max(0, int(u)));
      }
      const size_t g = c[0] + (size_t)gd[0] * (c[1] + (size_t)gd[1] * c[2]);
      for (uint32_t q = gstart[g]; q < gstart[g + 1]; q++) {
        const int e = gitems[q];
        V3 V[8]; float S[8];
        for (int i = 0; i < nv[e]; i++) { const int v = idx[ofs[e] + i]; V[i] = vtx[v]; S[i] = scal[v]; }
        float val;
        if (eltValue(nv[e], V, S, P, val)) return val;
      }
      return NAN;
    }
  };

  inline float fastLog(float f)
  {
    f = (f - 1.f) / (f + 1.f);
    const float f2 = f * f;
    float s = f;
    f *= f2; s += (1.f / 3.f) * f;
    f *= f2; s += (1.f / 5.f) * f;
    return s + s;
  }

  struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed * 0x9E3779B97F4A7C15ull + 1) {}
    float operator()() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return (s >> 40) * (1.f / 16777216.f); }
  };

  bool globalMaj = false, exactLog = false;

  // MCVolumeAccel::isProg for one part: returns hit t (or +inf)
  float tracePart(const Part &p, const XF &xf, V3 org, V3 dir, float tMax, Rng &rng, float *rgb)
  {
    // boxTest
    float t0 = 0.f, t1 = tMax;
    for (int d = 0; d < 3; d++) {
      const float inv = 1.f / get(dir, d);
      float a = (get(p.lo, d) - get(org, d)) * inv, b = (get(p.hi, d) - get(org, d)) * inv;
      if (a > b) std::swap(a, b);
      t0 = std::max(t0, a); t1 = std::min(t1, b);
    }
    if (t0 > t1) return INFINITY;
    float tUpper = t1;
    auto woodcock = [&](float lo_, float hi_, float maj) -> float {
      float t = lo_;
      while (true) {
        const float r = rng();
        const float dt = (exactLog ? -logf(1.f - r) : -fastLog(1.f - r)) / maj;
        t += dt;
        if (t >= hi_) return INFINITY;
        const V3 P = org + dir * t;
        const float s = p.sample(P);
        float c[3], w = 0.f;
        if (!std::isnan(s)) xf.map(s, c, w);
        if (w >= rng() * maj) { if (rgb) { rgb[0] = c[0]; rgb[1] = c[1]; rgb[2] = c[2]; } return t; }
      }
    };
    if (globalMaj) {
      if (p.globalMajorant == 0.f) return INFINITY;
      return woodcock(t0, tUpper, p.globalMajorant);
    }
    // dda3 in MC grid space (same arithmetic as DDA.h)
    V3 o = { (org.x - p.origin.x) / p.spacing.x, (org.y - p.origin.y) / p.spacing.y, (org.z - p.origin.z) / p.spacing.z };
    V3 dr = { dir.x / p.spacing.x, dir.y / p.spacing.y, dir.z / p.spacing.z };
    const V3 fs = { (float)p.dims[0], (float)p.dims[1], (float)p.dims[2] };
    V3 tnr, tfr;
    for (int d = 0; d < 3; d++) {
      const float a = (0.f - get(o, d)) / get(dr, d), b = (get(fs, d) - get(o, d)) / get(dr, d);
      set(tnr, d, std::min(a, b)); set(tfr, d, std::max(a, b));
    }
    const float ray_t0 = std::max(0.f, std::max(tnr.x, std::max(tnr.y, tnr.z)));
    const float ray_t1 = std::min(tUpper, std::min(tfr.x, std::min(tfr.y, tfr.z)));
    if (ray_t0 > ray_t1) return INFINITY;
    const V3 oin = o + dr * ray_t0;
    int cell[3], stop[3], delta[3];
    V3 tnext, tstep;
    for (int d = 0; d < 3; d++) {
      const float fc = std::max(0.f, std::min(get(fs, d) - 1.f, floorf(get(oin, d))));
      const float fe = get(dr, d) > 0.f ? fc + 1.f : fc;
      set(tstep, d, std::fabs(1.f / get(dr, d)));
      set(tnext, d, std::fabs(fe - get(oin, d)) * get(tstep, d));
      cell[d] = (int)fc;
      stop[d] = get(dr, d) > 0.f ? p.dims[d] : -1;
      delta[d] = get(dr, d) > 0.f ? 1 : -1;
    }
    float next_begin = 0.f;
    while (true) {
      const float tc = std::min(tnext.x, std::min(tnext.y, tnext.z));
      const float c0 = ray_t0 + next_begin, c1 = ray_t0 + std::min(tc, tUpper);
      const float maj = p.majorant[cell[0] + (size_t)p.dims[0] * (cell[1] + (size_t)p.dims[1] * cell[2])];
      if (maj != 0.f) {
        const float t = woodcock(c0, std::min(c1, tUpper), maj);
        if (t < INFINITY) return t;
      }
      next_begin = tc;
      for (int d = 0; d < 3; d++)
        if (get(tnext, d) == tc) {
          set(tnext, d, get(tnext, d) + get(tstep, d));
          cell[d] += delta[d];
          if (cell[d] == stop[d]) return INFINITY;
        }
    }
  }

} // namespace

int main(int ac, char **av)
{
  try {
    std::string xfFile, out;
    std::vector<std::string> files;
    float cam[9] = { 105, -105, 65, 0, 0, 0, 0, 0, 1 }, fovy = 40.f;
    int W = 960, H = 540, spp = 16;
    for (int i = 1; i < ac; i++) {
      const std::string a = av[i];
      if (a == "-xf") xfFile = av[++i];
      else if (a == "-o") out = av[++i];
      else if (a == "--camera") for (int q = 0; q < 9; q++) cam[q] = std::stof(av[++i]);
      else if (a == "-fovy") fovy = std::stof(av[++i]);
      else if (a == "-res") { W = std::stoi(av[++i]); H = std::stoi(av[++i]); }
      else if (a == "-spp") spp = std::stoi(av[++i]);
      else if (a == "--global-majorant") globalMaj = true;
      else if (a == "--exact-log") exactLog = true;
      else files.push_back(a);
    }
    if (xfFile.empty() || out.empty() || files.empty()) {
      std::cerr << "usage: barneyVolumeEmu -xf tf.xf -o prefix [--camera 9 floats] [-fovy f] [-res W H] [-spp n] "
                   "[--global-majorant] [--exact-log] parts.umesh...\n";
      return 1;
    }
    XF xf;
    xf.load(xfFile);
    std::vector<Part> parts(files.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < (int)files.size(); i++) parts[i].load(files[i], xf);
    size_t mcNonZero = 0, mcTotal = 0;
    for (auto &p : parts) { mcTotal += p.majorant.size(); for (float m : p.majorant) mcNonZero += m > 0.f; }
    printf("%zu parts loaded; MC grid e.g. %dx%dx%d; %zu of %zu MCs have majorant > 0\n", parts.size(),
           parts[0].dims[0], parts[0].dims[1], parts[0].dims[2], mcNonZero, mcTotal);

    V3 eye = { cam[0], cam[1], cam[2] }, at = { cam[3], cam[4], cam[5] }, up = { cam[6], cam[7], cam[8] };
    V3 dir = at - eye; dir = dir * (1.f / std::sqrt(dot(dir, dir)));
    V3 right = cross(dir, up); right = right * (1.f / std::sqrt(dot(right, right)));
    V3 cu = cross(right, dir);
    const float th = tanf(0.5f * fovy * (float)M_PI / 180.f), aspect = (float)W / H;
    std::vector<float> img((size_t)W * H * 3, 0.f), depth((size_t)W * H, 0.f);
#pragma omp parallel for schedule(dynamic, 2)
    for (int py = 0; py < H; py++)
      for (int px = 0; px < W; px++) {
        Rng rng(((uint64_t)py << 32) ^ (uint64_t)px ^ 0x1234567ull);
        float acc[3] = { 0, 0, 0 }, dacc = 0; int nh = 0;
        for (int s = 0; s < spp; s++) {
          const float sx = (2.f * (px + rng()) / W - 1.f) * th * aspect, sy = (1.f - 2.f * (py + rng()) / H) * th;
          V3 d = dir + right * sx + cu * sy;
          d = d * (1.f / std::sqrt(dot(d, d)));
          float tBest = INFINITY, rgbBest[3] = { 0, 0, 0 };
          for (auto &p : parts) {
            float rgb[3];
            const float t = tracePart(p, xf, eye, d, tBest, rng, rgb);
            if (t < tBest) { tBest = t; memcpy(rgbBest, rgb, sizeof(rgb)); }
          }
          if (tBest < INFINITY) { for (int q = 0; q < 3; q++) acc[q] += rgbBest[q]; dacc += tBest; nh++; }
        }
        for (int q = 0; q < 3; q++) img[((size_t)py * W + px) * 3 + q] = acc[q] / spp;
        depth[(size_t)py * W + px] = nh ? dacc / nh : 0.f;
      }
    // write PPMs: colour, and depth (near = bright) normalised to the hit range
    float dmin = 1e30f, dmax = 0;
    for (float d : depth) if (d > 0) { dmin = std::min(dmin, d); dmax = std::max(dmax, d); }
    FILE *f = fopen((out + "_color.ppm").c_str(), "wb");
    fprintf(f, "P6\n%d %d\n255\n", W, H);
    for (float v : img) fputc((int)std::min(255.f, 255.f * v), f);
    fclose(f);
    f = fopen((out + "_depth.ppm").c_str(), "wb");
    fprintf(f, "P6\n%d %d\n255\n", W, H);
    for (float d : depth) { const int g = d > 0 ? (int)(255.f * (1.f - 0.85f * (d - dmin) / (dmax - dmin))) : 0; fputc(g, f); fputc(g, f); fputc(g, f); }
    fclose(f);
    printf("wrote %s_{color,depth}.ppm (%dx%d, %d spp, %s, %s), hit depth [%g, %g]\n", out.c_str(), W, H, spp,
           globalMaj ? "one majorant per part" : "MC-grid majorants", exactLog ? "logf" : "fastLog", dmin, dmax);
  } catch (const std::exception &e) {
    std::cerr << "barneyVolumeEmu: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
