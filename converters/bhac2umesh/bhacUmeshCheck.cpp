// bhacUmeshCheck: validate .umesh files written by bhac2umesh.
//
// For every file: loads it with umesh::UMesh::loadFrom(), checks index ranges,
// NaN/Inf in vertices and per-vertex scalars, the element orientation with the
// criteria of umesh/apps/fixNegativeVolumeElements.cpp, the signed volumes, and
// that barney's point-in-element test finds points that lie inside the
// spherical cell. The Newton inversion below (hex and prism shape functions,
// VTK node order) is a host copy of
// barney/native/umesh/common/ElementIntersection.h (intersectHexEXT /
// intersectPrismEXT, adapted from OpenVKL, Apache-2.0).
//
// usage: bhacUmeshCheck [--expect r|z|torus] <file.umesh> [more.umesh ...]
//   --expect: the files were written with bhac2umesh --analytic <f>; compares the
//   value barney would interpolate at the sample points with the exact field
//   --bad-order-control: negative control for --expect, interpolates with the
//   scalars of hex corners 2 and 3 (wedge corners 0 and 1) swapped (must give large errors)

#include "umesh/UMesh.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using umesh::UMesh;
using umesh::vec3f;

namespace {

  struct V3 { double x, y, z; };
  inline V3 operator-(V3 a, V3 b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
  inline V3 operator+(V3 a, V3 b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
  inline V3 operator*(V3 a, double s) { return { a.x * s, a.y * s, a.z * s }; }
  inline double dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
  inline V3 cross(V3 a, V3 b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }
  inline double det(V3 a, V3 b, V3 c) { return dot(a, cross(b, c)); }
  inline V3 toV3(const vec3f &v) { return { v.x, v.y, v.z }; }

  double tetVol6(V3 v0, V3 v1, V3 v2, V3 v3) { return dot(v3 - v0, cross(v1 - v0, v2 - v0)); }

  // ---- barney/OpenVKL shape functions (VTK node order) ----
  void hexW(const double *p, double *sf, double *d)
  {
    const double rm = 1 - p[0], sm = 1 - p[1], tm = 1 - p[2];
    sf[0] = rm*sm*tm; sf[1] = p[0]*sm*tm; sf[2] = p[0]*p[1]*tm; sf[3] = rm*p[1]*tm;
    sf[4] = rm*sm*p[2]; sf[5] = p[0]*sm*p[2]; sf[6] = p[0]*p[1]*p[2]; sf[7] = rm*p[1]*p[2];
    d[0] = -sm*tm; d[1] = sm*tm; d[2] = p[1]*tm; d[3] = -p[1]*tm;
    d[4] = -sm*p[2]; d[5] = sm*p[2]; d[6] = p[1]*p[2]; d[7] = -p[1]*p[2];
    d[8] = -rm*tm; d[9] = -p[0]*tm; d[10] = p[0]*tm; d[11] = rm*tm;
    d[12] = -rm*p[2]; d[13] = -p[0]*p[2]; d[14] = p[0]*p[2]; d[15] = rm*p[2];
    d[16] = -rm*sm; d[17] = -p[0]*sm; d[18] = -p[0]*p[1]; d[19] = -rm*p[1];
    d[20] = rm*sm; d[21] = p[0]*sm; d[22] = p[0]*p[1]; d[23] = rm*p[1];
  }
  void prismW(const double *p, double *sf, double *d)
  {
    sf[0] = (1 - p[0] - p[1]) * (1 - p[2]); sf[1] = p[0] * (1 - p[2]); sf[2] = p[1] * (1 - p[2]);
    sf[3] = (1 - p[0] - p[1]) * p[2]; sf[4] = p[0] * p[2]; sf[5] = p[1] * p[2];
    d[0] = -1 + p[2]; d[1] = 1 - p[2]; d[2] = 0; d[3] = -p[2]; d[4] = p[2]; d[5] = 0;
    d[6] = -1 + p[2]; d[7] = 0; d[8] = 1 - p[2]; d[9] = -p[2]; d[10] = 0; d[11] = p[2];
    d[12] = -1 + p[0] + p[1]; d[13] = -p[0]; d[14] = -p[1]; d[15] = 1 - p[0] - p[1]; d[16] = p[0]; d[17] = p[1];
  }

  // Newton inversion as in barney; returns true if P is inside, pc = local coords
  bool locate(int n, const V3 *V, V3 P, double *pc, double *jacDet, double *w = nullptr)
  {
    double p[3] = { 0.5, 0.5, 0.5 }, sf[8], d[24];
    bool converged = false;
    for (int it = 0; it < 10; it++) {
      if (n == 8) hexW(p, sf, d); else prismW(p, sf, d);
      V3 f = { 0, 0, 0 }, r = { 0, 0, 0 }, s = { 0, 0, 0 }, t = { 0, 0, 0 };
      for (int i = 0; i < n; i++) {
        f = f + V[i] * sf[i]; r = r + V[i] * d[i]; s = s + V[i] * d[i + n]; t = t + V[i] * d[i + 2 * n];
      }
      f = f - P;
      const double dd = det(r, s, t);
      if (it == 0) *jacDet = dd;
      if (std::fabs(dd) < 1e-12) return false;
      const double d0 = det(f, s, t) / dd, d1 = det(r, f, t) / dd, d2 = det(r, s, f) / dd;
      p[0] -= d0; p[1] -= d1; p[2] -= d2;
      if (std::fabs(d0) < 1e-4 && std::fabs(d1) < 1e-4 && std::fabs(d2) < 1e-4) { converged = true; break; }
      if (std::fabs(p[0]) > 1e6 || std::fabs(p[1]) > 1e6 || std::fabs(p[2]) > 1e6) return false;
    }
    pc[0] = p[0]; pc[1] = p[1]; pc[2] = p[2];
    if (!converged) return false;
    if (w) { if (n == 8) hexW(p, sf, d); else prismW(p, sf, d); for (int i = 0; i < n; i++) w[i] = sf[i]; }
    const double lo = -1e-6, hi = 1 + 1e-6;
    bool in = p[0] >= lo && p[0] <= hi && p[1] >= lo && p[1] <= hi && p[2] >= lo && p[2] <= hi;
    if (n == 6) in = in && p[0] + p[1] <= hi;
    return in;
  }

  // a point inside the spherical cell spanned by the element's vertices: the
  // (r, theta, phi) mid point (+ offsets f in [0,1]) of the vertex ranges
  V3 sphericalInterior(int n, const V3 *V, const double f[3])
  {
    double rmin = 1e30, rmax = 0, tmin = 10, tmax = -10, pmin = 1e30, pmax = -1e30, pref = 0;
    bool haveRef = false;
    for (int i = 0; i < n; i++) {
      const double r = std::sqrt(dot(V[i], V[i]));
      const double t = std::acos(std::max(-1.0, std::min(1.0, V[i].z / r)));
      rmin = std::min(rmin, r); rmax = std::max(rmax, r);
      tmin = std::min(tmin, t); tmax = std::max(tmax, t);
      const double rho = std::sqrt(V[i].x * V[i].x + V[i].y * V[i].y);
      if (rho < 1e-6 * r) continue;   // on the axis: phi undefined
      double p = std::atan2(V[i].y, V[i].x);
      if (!haveRef) { pref = p; haveRef = true; }
      while (p - pref > M_PI) p -= 2 * M_PI;
      while (p - pref < -M_PI) p += 2 * M_PI;
      pmin = std::min(pmin, p); pmax = std::max(pmax, p);
    }
    const double r = rmin + f[0] * (rmax - rmin), t = tmin + f[1] * (tmax - tmin), p = pmin + f[2] * (pmax - pmin);
    return { r * std::sin(t) * std::cos(p), r * std::sin(t) * std::sin(p), r * std::cos(t) };
  }

  // analytic test fields of bhac2umesh --analytic
  std::string expect;
  double analytic(V3 P)
  {
    const double r = std::sqrt(dot(P, P)), R = std::sqrt(P.x * P.x + P.y * P.y);
    if (expect == "r") return r;
    if (expect == "z") return P.z;
    const double d2 = (R - 20.0) * (R - 20.0) + P.z * P.z;
    return std::max(-9.0, -d2 / (2.0 * 16.0) / std::log(10.0));
  }
  double maxErr = 0, maxRelErr = 0, sumErr = 0, sumSpan = 0;
  size_t numSamples = 0;
  bool badOrderControl = false;
  std::string worst;   // negative control: swap the scalars of hex corners 2 and 3

  // value barney would return at P (weights at the local coords) vs the analytic field
  void compareValue(int n, const int *idx, const V3 *V, V3 P, const UMesh &m)
  {
    double pc[3], jd, w[8];
    if (!locate(n, V, P, pc, &jd, w)) return;
    double val = 0, vmin = 1e30, vmax = -1e30;
    for (int i = 0; i < n; i++) {
      int ii = i;
      if (badOrderControl && n == 8 && (i == 2 || i == 3)) ii = 5 - i;
      if (badOrderControl && n == 6 && (i == 0 || i == 1)) ii = 1 - i;
      const double s = m.perVertex->values[idx[ii]];
      val += w[i] * s; vmin = std::min(vmin, s); vmax = std::max(vmax, s);
    }
    const double err = std::fabs(val - analytic(P));
    if (err > maxErr) {
      char buf[512];
      snprintf(buf, sizeof(buf), "n=%d P=(%g %g %g) |P|=%g pc=(%g %g %g) interp=%g exact=%g vals[%g..%g]",
               n, P.x, P.y, P.z, std::sqrt(dot(P, P)), pc[0], pc[1], pc[2], val, analytic(P), vmin, vmax);
      worst = buf;
    }
    maxErr = std::max(maxErr, err);
    sumErr += err; sumSpan += vmax - vmin; numSamples++;
    if (vmax - vmin > 1e-6 * std::max(1.0, std::fabs(vmax))) maxRelErr = std::max(maxRelErr, err / (vmax - vmin));
  }

  bool checkFile(const std::string &fn)
  {
    UMesh::SP m = UMesh::loadFrom(fn);
    const size_t nv = m->vertices.size();
    size_t errors = 0;
    if (!m->perVertex || m->perVertex->values.size() != nv) {
      printf("%s: per-vertex scalar missing or of wrong size\n", fn.c_str());
      return false;
    }
    size_t nanV = 0, nanS = 0;
    for (auto &v : m->vertices) if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z)) nanV++;
    for (float s : m->perVertex->values) if (!std::isfinite(s)) nanS++;
    auto badIdx = [&](int i) { return i < 0 || (size_t)i >= nv; };
    size_t badIndex = 0, badOrient = 0, negVol = 0, notFound = 0, samples = 0;
    double volume = 0, minVol = 1e30, maxVol = 0, minJac = 1e30, maxJac = -1e30;
    std::vector<char> used(nv, 0);
    const double fr[3][3] = { { 0.5, 0.5, 0.5 }, { 0.3, 0.7, 0.4 }, { 0.8, 0.25, 0.65 } };

    for (auto &h : m->hexes) {
      V3 V[8];
      bool ok = true;
      for (int i = 0; i < 8; i++) { if (badIdx(h[i])) ok = false; else { V[i] = toV3(m->vertices[h[i]]); used[h[i]] = 1; } }
      if (!ok) { badIndex++; continue; }
      // umesh orientation criterion (fixNegativeVolumeElements)
      V3 c = { 0, 0, 0 }, b = { 0, 0, 0 };
      for (int i = 0; i < 8; i++) c = c + V[i] * 0.125;
      for (int i = 0; i < 4; i++) b = b + V[i] * 0.25;
      if (tetVol6(V[0], V[1], b, c) < 0) badOrient++;
      static const int t[6][4] = { {0,1,2,6}, {0,2,3,6}, {0,3,7,6}, {0,7,4,6}, {0,4,5,6}, {0,5,1,6} };
      double vol = 0;
      for (auto &q : t) vol += tetVol6(V[q[0]], V[q[1]], V[q[2]], V[q[3]]) / 6.0;
      if (vol <= 0) negVol++;
      volume += vol; minVol = std::min(minVol, vol); maxVol = std::max(maxVol, vol);
      for (auto &f : fr) {
        double pc[3], jd;
        samples++;
        const V3 P = sphericalInterior(8, V, f);
        if (!locate(8, V, P, pc, &jd)) notFound++;
        minJac = std::min(minJac, jd); maxJac = std::max(maxJac, jd);
        if (!expect.empty()) { int idx[8]; for (int i = 0; i < 8; i++) idx[i] = h[i]; compareValue(8, idx, V, P, *m); }
      }
    }
    double wMinJac = 1e30, wMaxJac = -1e30;
    for (auto &w : m->wedges) {
      V3 V[6];
      bool ok = true;
      for (int i = 0; i < 6; i++) { if (badIdx(w[i])) ok = false; else { V[i] = toV3(m->vertices[w[i]]); used[w[i]] = 1; } }
      if (!ok) { badIndex++; continue; }
      V3 b = (V[0] + V[1] + V[3] + V[4]) * 0.25;
      if (tetVol6(V[3], V[4], V[5], b) < 0) badOrient++;
      // prism = tets (3,4,5,0), (3,1,4,0)... computed as the degenerate hex (3,4,5,5 | 0,1,2,2)
      const V3 H[8] = { V[3], V[4], V[5], V[5], V[0], V[1], V[2], V[2] };
      static const int t[6][4] = { {0,1,2,6}, {0,2,3,6}, {0,3,7,6}, {0,7,4,6}, {0,4,5,6}, {0,5,1,6} };
      double vol = 0;
      for (auto &q : t) vol += tetVol6(H[q[0]], H[q[1]], H[q[2]], H[q[3]]) / 6.0;
      if (vol <= 0) negVol++;
      volume += vol; minVol = std::min(minVol, vol); maxVol = std::max(maxVol, vol);
      for (auto &f : fr) {
        double pc[3], jd;
        samples++;
        const V3 P = sphericalInterior(6, V, f);
        if (!locate(6, V, P, pc, &jd)) notFound++;
        wMinJac = std::min(wMinJac, jd); wMaxJac = std::max(wMaxJac, jd);
        if (!expect.empty()) { int idx[6]; for (int i = 0; i < 6; i++) idx[i] = w[i]; compareValue(6, idx, V, P, *m); }
      }
    }
    size_t badTet = 0;
    for (auto &t : m->tets) {
      bool ok = true;
      V3 V[4];
      for (int i = 0; i < 4; i++) { if (badIdx(t[i])) ok = false; else { V[i] = toV3(m->vertices[t[i]]); used[t[i]] = 1; } }
      if (!ok) { badIndex++; continue; }
      // barney tetScalar: the centroid must be on the positive side of all 4 planes
      const double vol = tetVol6(V[0], V[1], V[2], V[3]) / 6.0;
      if (vol <= 0) { negVol++; badTet++; }
      volume += vol; minVol = std::min(minVol, vol); maxVol = std::max(maxVol, vol);
    }
    size_t unused = 0;
    for (char u : used) if (!u) unused++;
    const size_t other = m->pyrs.size() + m->grids.size() + m->polyOffsets.size() +
                         m->triangles.size() + m->quads.size();
    errors = nanV + nanS + badIndex + badOrient + negVol + notFound + other;
    auto bb = m->getBounds();
    printf("%s: %zu tets (%zu non-positive), %zu hexes, %zu wedges, %zu other prims, %zu vertices (%zu unused), bounds (%g %g %g)-(%g %g %g), "
           "values [%g, %g]\n  NaN vertices %zu, NaN scalars %zu, bad indices %zu, bad orientation %zu, "
           "non-positive volume %zu, volume %.6g (min %.3g max %.3g)\n"
           "  barney point location: %zu of %zu interior samples not found; initial Jacobian det hex [%.3g, %.3g], "
           "wedge [%.3g, %.3g] -> %s\n",
           fn.c_str(), m->tets.size(), badTet, m->hexes.size(), m->wedges.size(), other, nv, unused,
           bb.lower.x, bb.lower.y, bb.lower.z, bb.upper.x, bb.upper.y, bb.upper.z,
           m->perVertex->valueRange.lower, m->perVertex->valueRange.upper,
           nanV, nanS, badIndex, badOrient, negVol, volume, minVol, maxVol, notFound, samples,
           minJac, maxJac, wMinJac, wMaxJac, errors ? "FAILED" : "OK");
    if (!expect.empty()) {
      printf("  analytic '%s': max |interpolated - exact| %.4g, max error / element value span %.4g\n",
             expect.c_str(), maxErr, maxRelErr);
      printf("  mean |error| %.4g, mean element value span %.4g, ratio %.4g\n", sumErr / std::max<size_t>(1, numSamples),
             sumSpan / std::max<size_t>(1, numSamples), sumErr / std::max(1e-300, sumSpan));
      printf("  worst sample: %s\n", worst.c_str());
      sumErr = sumSpan = 0; numSamples = 0;
      maxErr = maxRelErr = 0;
    }
    return errors == 0;
  }

} // namespace

int main(int ac, char **av)
{
  if (ac < 2) {
    std::cerr << "usage: bhacUmeshCheck <file.umesh> [more.umesh ...]" << std::endl;
    return 2;
  }
  int failed = 0, numFiles = 0;
  for (int i = 1; i < ac; i++) {
    if (std::string(av[i]) == "--expect" && i + 1 < ac) { expect = av[++i]; continue; }
    if (std::string(av[i]) == "--bad-order-control") { badOrderControl = true; continue; }
    numFiles++;
    try {
      if (!checkFile(av[i])) failed++;
    } catch (const std::exception &e) {
      std::cerr << av[i] << ": " << e.what() << std::endl;
      failed++;
    }
  }
  printf("%d of %d files failed\n", failed, numFiles);
  return failed ? 1 : 0;
}
