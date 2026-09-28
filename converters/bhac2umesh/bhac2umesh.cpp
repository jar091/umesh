// bhac2umesh: convert a BHAC (MPI-AMRVAC based GRMHD code) snapshot on a
// spherical modified-Kerr-Schild (MKS) grid into unstructured meshes (.umesh)
// that haystack/barney can render, split into N spatially coherent parts for
// data-parallel rendering.
//
//  - every cell becomes a hexahedron whose 8 corners are the cell-corner code
//    coordinates (x1, x2, x3) mapped MKS -> (r, theta, phi) -> Cartesian;
//    cells touching the polar axis (theta = 0 or pi) become wedges, because two
//    of their corner pairs collapse onto the axis
//  - vertices are shared between neighbouring cells, across the phi seam
//    (phi = 0 == 2 pi) and on the polar axis (one vertex per radius and pole)
//  - one per-vertex scalar: the average of the (transformed) values of all
//    cells sharing the vertex
//  - element vertex order is the umesh/VTK order (see writeCell() below)
//
// Only single-level snapshots (levmax = 1, i.e. mxnest = 1) are supported: AMR
// levels would need hanging-node handling.

#include "BhacSnapshot.h"
#include "umesh/UMesh.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using umesh::UMesh;
using umesh::vec3f;

namespace {

  // --------------------------------------------------------------------------
  // command line
  // --------------------------------------------------------------------------
  struct Options {
    std::string datFile, parFile, outBase;
    int numParts = 1;
    std::string var = "rho";
    std::string scale = "auto";   // auto | log | lin
    double rmin = 0.0, rmax = 0.0;
    double mksH = 0.35, mksR0 = 0.0;
    double spin = NAN, mass = NAN, gammaAd = NAN;   // NAN: from eqpar of the snapshot
    bool listVars = false;
    std::string analytic;         // test field instead of the snapshot data: r | z | torus
    bool tets = false;            // split every cell into tets (6 per hex, 3 per axis cell)
    // CPU reference render (diagnostics, no .umesh output)
    std::string refRender;
    double refIso = -3.4, refFovy = 40.0, refRmax = 200.0, refStep = 0.05;
    double refCam[9] = { 105, -105, 65, 0, 0, 0, 0, 0, 1 };
    int refRes[2] = { 960, 540 };
  };

  void usage(const std::string &err = "")
  {
    if (!err.empty()) std::cerr << "Error: " << err << "\n\n";
    std::cout <<
      "usage: bhac2umesh <dataNNNN.dat> --par <amrvac.par> -o <outBase> [options]\n"
      "  writes <outBase>_<i>.umesh (i = 0..n-1) and <outBase>.json\n"
      "options:\n"
      "  -n <parts>        number of spatially coherent parts (default 1)\n"
      "  --var <name>      scalar: rho (default, = d/lfac), p, temp (= p/rho), bsq (= b^2),\n"
      "                    sigma (= b^2/rho), beta (= 2p/b^2), or a stored variable\n"
      "                    (d s1 s2 s3 tau b1 b2 b3 Ds dtr1 lfac xi)\n"
      "  --log | --lin     log10 or linear scale (default: log10 for rho p temp bsq sigma beta,\n"
      "                    linear otherwise); non-positive values get the smallest positive value\n"
      "  --rmin <r>        drop radial shells whose cell centre has r < rmin (default: keep all)\n"
      "  --rmax <r>        drop radial shells whose cell centre has r > rmax (default: keep all)\n"
      "  --mks-h <h>       MKS theta squeeze coordpar(h_) (default 0.35; not stored in the .dat)\n"
      "  --mks-r0 <R0>     MKS radial offset coordpar(R0_) (default 0)\n"
      "  --spin <a> --mass <M> --gamma <g>   override eqpar(a_), eqpar(m_), eqpar(gamma_)\n"
      "  --list-vars       print the snapshot header and the variables, then exit\n"
      "  --tets            write tetrahedra instead of hexes/wedges: every cell split into the 6\n"
      "                    tets around its 0-6 diagonal (conforming between neighbours, positive\n"
      "                    orientation); axis cells give 3 tets\n"
      "  --analytic <f>    test mode: same grid/geometry/averaging, but the cell values are an\n"
      "                    analytic function of the Cartesian cell centre instead of the data:\n"
      "                    r (= |x|), z, torus (= log10 of a Gaussian torus,\n"
      "                    -((R_cyl-20)^2+z^2)/(2*4^2)/ln10, floored at -9; linear scale)\n"
      "  --ref-render <out.ppm>  diagnostics: instead of writing .umesh, ray-cast the iso-surface\n"
      "                    --ref-iso <v> (default -3.4) of the per-vertex field, interpolated\n"
      "                    trilinearly in the logical cell (what the hexes represent), Lambert\n"
      "                    shaded; --ref-camera px py pz lx ly lz ux uy uz, --ref-fovy, --ref-res W H,\n"
      "                    --ref-rmax (march inside r < rmax), --ref-step (march step, M)\n";
    exit(err.empty() ? 0 : 1);
  }

  Options parseArgs(int ac, char **av)
  {
    Options o;
    for (int i = 1; i < ac; i++) {
      const std::string a = av[i];
      auto next = [&]() -> std::string {
        if (i + 1 >= ac) usage("missing value after " + a);
        return av[++i];
      };
      if (a == "-h" || a == "--help") usage();
      else if (a == "--par") o.parFile = next();
      else if (a == "-o") o.outBase = next();
      else if (a == "-n") o.numParts = std::stoi(next());
      else if (a == "--var") o.var = next();
      else if (a == "--log") o.scale = "log";
      else if (a == "--lin") o.scale = "lin";
      else if (a == "--rmin") o.rmin = std::stod(next());
      else if (a == "--rmax") o.rmax = std::stod(next());
      else if (a == "--mks-h") o.mksH = std::stod(next());
      else if (a == "--mks-r0") o.mksR0 = std::stod(next());
      else if (a == "--spin") o.spin = std::stod(next());
      else if (a == "--mass") o.mass = std::stod(next());
      else if (a == "--gamma") o.gammaAd = std::stod(next());
      else if (a == "--list-vars") o.listVars = true;
      else if (a == "--analytic") o.analytic = next();
      else if (a == "--tets") o.tets = true;
      else if (a == "--ref-render") o.refRender = next();
      else if (a == "--ref-iso") o.refIso = std::stod(next());
      else if (a == "--ref-fovy") o.refFovy = std::stod(next());
      else if (a == "--ref-rmax") o.refRmax = std::stod(next());
      else if (a == "--ref-step") o.refStep = std::stod(next());
      else if (a == "--ref-res") { o.refRes[0] = std::stoi(next()); o.refRes[1] = std::stoi(next()); }
      else if (a == "--ref-camera") { for (int q = 0; q < 9; q++) o.refCam[q] = std::stod(next()); }
      else if (!a.empty() && a[0] == '-') usage("unknown option " + a);
      else if (o.datFile.empty()) o.datFile = a;
      else usage("more than one input file given");
    }
    if (o.datFile.empty()) usage("no input .dat file");
    if (o.parFile.empty()) usage("--par is required (the grid is not stored in the .dat file)");
    if (!o.listVars && o.refRender.empty() && o.outBase.empty()) usage("no output base (-o)");
    if (o.numParts < 1) usage("-n must be >= 1");
    return o;
  }

  // --------------------------------------------------------------------------
  // geometry helpers
  // --------------------------------------------------------------------------
  struct Box { int lo[3], hi[3]; size_t numCells() const {
      return (size_t)(hi[0]-lo[0])*(size_t)(hi[1]-lo[1])*(size_t)(hi[2]-lo[2]); } };

  inline double dot3(const double *a, const double *b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
  inline void sub3(const vec3f &a, const vec3f &b, double *r) { r[0]=(double)a.x-b.x; r[1]=(double)a.y-b.y; r[2]=(double)a.z-b.z; }
  inline void cross3(const double *a, const double *b, double *r)
  { r[0]=a[1]*b[2]-a[2]*b[1]; r[1]=a[2]*b[0]-a[0]*b[2]; r[2]=a[0]*b[1]-a[1]*b[0]; }

  // umesh's tet "volume": dot(v3-v0, cross(v1-v0, v2-v0)) (6x signed volume),
  // as in umesh/apps/fixNegativeVolumeElements.cpp
  double tetVol6(const vec3f &v0, const vec3f &v1, const vec3f &v2, const vec3f &v3)
  {
    double a[3], b[3], c[3], n[3];
    sub3(v1, v0, a); sub3(v2, v0, b); sub3(v3, v0, c);
    cross3(a, b, n);
    return dot3(c, n);
  }
  vec3f avg(std::initializer_list<vec3f> l)
  {
    double s[3] = { 0, 0, 0 };
    for (auto &v : l) { s[0] += v.x; s[1] += v.y; s[2] += v.z; }
    const double n = (double)l.size();
    return vec3f((float)(s[0]/n), (float)(s[1]/n), (float)(s[2]/n));
  }

  // Orientation criteria of umesh/apps/fixNegativeVolumeElements.cpp (the
  // tool umesh provides to repair element orientation): true = correct
  bool hexOrientedOK(const vec3f *v)
  {
    const vec3f c = avg({ v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7] });
    const vec3f b = avg({ v[0], v[1], v[2], v[3] });
    return tetVol6(v[0], v[1], b, c) >= 0.0;
  }
  bool wedgeOrientedOK(const vec3f *v)
  {
    const vec3f b = avg({ v[0], v[1], v[3], v[4] });
    return tetVol6(v[3], v[4], v[5], b) >= 0.0;
  }
  // signed volume of a hex from 6 tets around the diagonal 0-6 (VTK order;
  // positive for the VTK/umesh orientation)
  double hexVolume(const vec3f *v)
  {
    static const int t[6][4] = { {0,1,2,6}, {0,2,3,6}, {0,3,7,6}, {0,7,4,6}, {0,4,5,6}, {0,5,1,6} };
    double s = 0.0;
    for (auto &q : t) s += tetVol6(v[q[0]], v[q[1]], v[q[2]], v[q[3]]);
    return s / 6.0;
  }
  // volume of a wedge in umesh/VTK order (front (0,1,2) faces away from (3,4,5));
  // computed on the equivalent degenerate hex (3,4,5,5 | 0,1,2,2), positive for
  // a correctly oriented wedge
  double wedgeVolume(const vec3f *v)
  {
    const vec3f h[8] = { v[3], v[4], v[5], v[5], v[0], v[1], v[2], v[2] };
    return hexVolume(h);
  }

  // --------------------------------------------------------------------------
  // the converter
  // --------------------------------------------------------------------------
  struct Converter {
    Options opt;
    bhac::Snapshot snap;

    double spin = 0.0, mass = 1.0, gammaAd = 5.0/3.0;
    int N[3];                               // global cells (x1, x2, x3)
    std::vector<double> x1e, x2e, x3e;      // cell-edge code coordinates
    std::vector<double> re, the, phe;       // cell-edge r, theta, phi
    bool northPole = false, southPole = false, periodicPhi = false;
    int i0 = 0, i1 = 0;                     // kept radial cell range [i0, i1)

    // vertex numbering
    int jlo = 0, jhi = 0, nJreg = 0, nKv = 0;
    size_t numRegular = 0, numVerts = 0;

    std::vector<double> cellVal;            // transformed cell values, index cid()
    std::vector<vec3f> vertPos;
    std::vector<float> vertVal;
    bool logScale = true;
    double logFloor = 0.0;
    size_t numNonPositive = 0, numNonFinite = 0;
    double rawMin = 0, rawMax = 0;

    size_t cid(int i, int j, int k) const { return (size_t)i + (size_t)N[0] * ((size_t)j + (size_t)N[1] * (size_t)k); }

    double thetaOf(double x2) const { return x2 + 0.5 * opt.mksH * std::sin(2.0 * x2); }
    double rOf(double x1) const { return opt.mksR0 + std::exp(x1); }

    /*! global vertex id of logical corner (i, j, k), i in [0,N0], j in [0,N1], k in [0,N2] */
    size_t vid(int i, int j, int k) const
    {
      if (northPole && j == 0) return numRegular + (size_t)i;
      if (southPole && j == N[1]) return numRegular + (northPole ? (size_t)(N[0] + 1) : 0) + (size_t)i;
      const int kk = periodicPhi ? (k % N[2]) : k;
      return ((size_t)i * nJreg + (size_t)(j - jlo)) * nKv + (size_t)kk;
    }

    void setup();
    void computeCellValues();
    void computeVertices();
    void partition(const Box &b, int n, std::vector<Box> &out) const;
    void run();
    double sampleLogical(double x, double y, double z) const;
    void refRender() const;
  };

  void Converter::setup()
  {
    snap.open(opt.datFile, opt.parFile);
    if (snap.ndim != 3 || snap.typeaxial != "spherical")
      throw std::runtime_error("only 3D snapshots with typeaxial = 'spherical' (MKS) are supported");
    for (auto &l : snap.leaves)
      if (l.level != 1)
        throw std::runtime_error("the snapshot has refined blocks (levmax = " + std::to_string(snap.levmax) +
                                 "); only single-level grids are supported");
    // eqpar: gamma_=1, adiab_=2, m_=3, a_=4 (amrvacpar.f of the GR MHD physics)
    gammaAd = !std::isnan(opt.gammaAd) ? opt.gammaAd : (snap.neqpar >= 1 ? snap.eqpar[0] : 5.0/3.0);
    mass    = !std::isnan(opt.mass)    ? opt.mass    : (snap.neqpar >= 3 ? snap.eqpar[2] : 1.0);
    spin    = !std::isnan(opt.spin)    ? opt.spin    : (snap.neqpar >= 4 ? snap.eqpar[3] : 0.0);

    for (int d = 0; d < 3; d++) N[d] = snap.nxlone[d];
    auto edges = [&](int d, std::vector<double> &e) {
      e.resize(N[d] + 1);
      const double dx = (snap.xprobmax[d] - snap.xprobmin[d]) / N[d];
      for (int i = 0; i <= N[d]; i++) e[i] = snap.xprobmin[d] + i * dx;
      e[N[d]] = snap.xprobmax[d];
    };
    edges(0, x1e); edges(1, x2e); edges(2, x3e);
    re.resize(x1e.size()); the.resize(x2e.size()); phe = x3e;
    for (size_t i = 0; i < x1e.size(); i++) re[i] = rOf(x1e[i]);
    for (size_t j = 0; j < x2e.size(); j++) the[j] = thetaOf(x2e[j]);

    const double eps = 1e-9;
    northPole = std::fabs(the[0]) < eps;
    southPole = std::fabs(the[N[1]] - M_PI) < eps;
    if (northPole) the[0] = 0.0;
    if (southPole) the[N[1]] = M_PI;
    periodicPhi = std::fabs((x3e[N[2]] - x3e[0]) - 2.0 * M_PI) < 1e-9;
    if (N[1] < 2 && northPole && southPole)
      throw std::runtime_error("need at least 2 cells in theta");

    // radial cut by cell-centre radius
    i0 = 0; i1 = N[0];
    for (int i = 0; i < N[0]; i++) {
      const double rc = rOf(0.5 * (x1e[i] + x1e[i + 1]));
      if (opt.rmin > 0.0 && rc < opt.rmin) i0 = i + 1;
      if (opt.rmax > 0.0 && rc > opt.rmax) { i1 = i; break; }
    }
    if (i1 <= i0) throw std::runtime_error("no cells left after the --rmin/--rmax cut");

    jlo = northPole ? 1 : 0;
    jhi = southPole ? N[1] - 1 : N[1];
    nJreg = jhi - jlo + 1;
    nKv = periodicPhi ? N[2] : N[2] + 1;
    numRegular = (size_t)(N[0] + 1) * nJreg * nKv;
    numVerts = numRegular + (northPole ? N[0] + 1 : 0) + (southPole ? N[0] + 1 : 0);

    printf("BHAC snapshot %s: it %d, t %g, nw %d, nws %d, block %dx%dx%d, %d leaf blocks (%dx%dx%d), levmax %d\n",
           opt.datFile.c_str(), snap.iteration, snap.time, snap.nw, snap.nws, snap.nx[0], snap.nx[1], snap.nx[2],
           snap.nleafs, snap.ng[0], snap.ng[1], snap.ng[2], snap.levmax);
    printf("variables:");
    for (auto &w : snap.wnames) printf(" %s", w.c_str());
    printf("\neqpar (%d):", snap.neqpar);
    for (double e : snap.eqpar) printf(" %g", e);
    printf("\ngrid %dx%dx%d cells, MKS h %g R0 %g: x1 [%g, %g] -> r [%g, %g], theta [%g, %g], phi [%g, %g]\n",
           N[0], N[1], N[2], opt.mksH, opt.mksR0, x1e[0], x1e[N[0]], re[0], re[N[0]], the[0], the[N[1]],
           phe[0], phe[N[2]]);
    printf("a = %g, M = %g, gamma = %g, horizon r+ = %g; poles: north %d south %d, periodic phi %d\n",
           spin, mass, gammaAd, mass + std::sqrt(std::max(0.0, mass*mass - spin*spin)),
           (int)northPole, (int)southPole, (int)periodicPhi);
    printf("kept radial cells [%d, %d) -> r [%g, %g]\n", i0, i1, re[i0], re[i1]);
  }

  void Converter::computeCellValues()
  {
    if (!opt.analytic.empty()) {
      const std::string a = bhac::lower(opt.analytic);
      if (a != "r" && a != "z" && a != "torus")
        throw std::runtime_error("unknown --analytic field '" + opt.analytic + "' (r, z, torus)");
      opt.var = "analytic_" + a;
      logScale = false;
      cellVal.assign((size_t)N[0] * N[1] * N[2], 0.0);
      for (int k = 0; k < N[2]; k++)
        for (int j = 0; j < N[1]; j++)
          for (int i = i0; i < i1; i++) {
            // Cartesian cell centre (logical centre mapped like the vertices)
            const double r = rOf(0.5 * (x1e[i] + x1e[i + 1]));
            const double th = thetaOf(0.5 * (x2e[j] + x2e[j + 1]));
            const double R = r * std::sin(th), z = r * std::cos(th);
            double v;
            if (a == "r") v = r;
            else if (a == "z") v = z;
            else {
              const double d2 = (R - 20.0) * (R - 20.0) + z * z;
              v = std::max(-9.0, -d2 / (2.0 * 16.0) / std::log(10.0));
            }
            cellVal[cid(i, j, k)] = v;
          }
      rawMin = INFINITY; rawMax = -INFINITY;
      for (int k = 0; k < N[2]; k++)
        for (int j = 0; j < N[1]; j++)
          for (int i = i0; i < i1; i++) {
            rawMin = std::min(rawMin, cellVal[cid(i, j, k)]);
            rawMax = std::max(rawMax, cellVal[cid(i, j, k)]);
          }
      printf("analytic test field %s: cell range [%g, %g]\n", a.c_str(), rawMin, rawMax);
      return;
    }
    const std::string var = bhac::lower(opt.var);
    enum { RAW, RHO, P, TEMP, BSQ, SIGMA, BETA } kind = RAW;
    int rawIdx = -1;
    if (var == "rho") kind = RHO;
    else if (var == "p") kind = P;
    else if (var == "temp") kind = TEMP;
    else if (var == "bsq") kind = BSQ;
    else if (var == "sigma") kind = SIGMA;
    else if (var == "beta") kind = BETA;
    else {
      rawIdx = snap.varIndex(opt.var);
      if (rawIdx < 0) throw std::runtime_error("unknown variable '" + opt.var + "'");
    }
    logScale = (opt.scale == "log") || (opt.scale == "auto" && kind != RAW);

    const int iD = snap.varIndex("d"), iLfac = snap.varIndex("lfac"), iXi = snap.varIndex("xi");
    const int iS[3] = { snap.varIndex("s1"), snap.varIndex("s2"), snap.varIndex("s3") };
    const int iB[3] = { snap.varIndex("b1"), snap.varIndex("b2"), snap.varIndex("b3") };
    if (kind != RAW && (iD < 0 || iLfac < 0))
      throw std::runtime_error("derived variables need d and lfac");
    if ((kind == P || kind == TEMP || kind == BETA) && iXi < 0)
      throw std::runtime_error("pressure needs xi");
    if ((kind == BSQ || kind == SIGMA || kind == BETA) && (iXi < 0 || iS[0] < 0 || iS[1] < 0 || iS[2] < 0 ||
                                                           iB[0] < 0 || iB[1] < 0 || iB[2] < 0))
      throw std::runtime_error("b^2 needs xi, s1..s3, b1..b3");

    cellVal.assign((size_t)N[0] * N[1] * N[2], 0.0);
    std::ifstream f(snap.datFile, std::ios::binary);
    std::vector<double> w;
    const size_t nc = snap.ncellBlock;
    const int *nx = snap.nx;
    for (size_t b = 0; b < snap.leaves.size(); b++) {
      const bhac::Leaf &leaf = snap.leaves[b];
      const int I0 = (leaf.ig[0] - 1) * nx[0], J0 = (leaf.ig[1] - 1) * nx[1], K0 = (leaf.ig[2] - 1) * nx[2];
      if (I0 + nx[0] <= i0 || I0 >= i1) continue;   // block outside the radial cut
      snap.readBlock(f, b, w);
      for (int k = 0; k < nx[2]; k++)
        for (int j = 0; j < nx[1]; j++)
          for (int i = 0; i < nx[0]; i++) {
            const size_t c = (size_t)i + (size_t)nx[0] * ((size_t)j + (size_t)nx[1] * (size_t)k);
            auto W = [&](int v) { return w[(size_t)v * nc + c]; };
            const int I = I0 + i, J = J0 + j, K = K0 + k;
            double val;
            if (kind == RAW) {
              val = W(rawIdx);
            } else {
              const double lfac = W(iLfac);
              const double rho = W(iD) / (lfac > 0.0 ? lfac : 1.0);   // D = Gamma rho (not densitized)
              double p = 0.0, bsq = 0.0;
              if (kind == P || kind == TEMP || kind == BETA)
                // xi = Gamma^2 rho h, rho h = rho + gamma/(gamma-1) p (eos gamma)
                p = (gammaAd - 1.0) / gammaAd * (W(iXi) / (lfac * lfac) - rho);
              if (kind == BSQ || kind == SIGMA || kind == BETA) {
                // b^2 = B^2 / Gamma^2 + (B^i v_i)^2 with B^i v_i = S_i B^i / xi,
                // B^2 = gamma_ij B^i B^j with the Kerr-Schild 3-metric in MKS
                // coordinates (mod_metric.f get_g_component)
                const double x1 = 0.5 * (x1e[I] + x1e[I + 1]), x2 = 0.5 * (x2e[J] + x2e[J + 1]);
                const double r = rOf(x1), th = thetaOf(x2), e1 = std::exp(x1);
                const double st = std::sin(th), ct = std::cos(th);
                const double sig = r * r + spin * spin * ct * ct;
                const double z = 2.0 * mass * r / sig;
                const double dth = 1.0 + opt.mksH * std::cos(2.0 * x2);
                const double g11 = e1 * e1 * (1.0 + z);
                const double g22 = sig * dth * dth;
                const double g33 = st * st * (spin * spin + r * r + 2.0 * spin * spin * mass * r * st * st / sig);
                const double g13 = -e1 * spin * (1.0 + z) * st * st;
                const double B1 = W(iB[0]), B2 = W(iB[1]), B3 = W(iB[2]);
                const double B2sq = g11 * B1 * B1 + g22 * B2 * B2 + g33 * B3 * B3 + 2.0 * g13 * B1 * B3;
                const double SB = W(iS[0]) * B1 + W(iS[1]) * B2 + W(iS[2]) * B3;
                const double vB = SB / W(iXi);
                bsq = B2sq / (lfac * lfac) + vB * vB;
              }
              switch (kind) {
              case RHO: val = rho; break;
              case P: val = p; break;
              case TEMP: val = p / rho; break;
              case BSQ: val = bsq; break;
              case SIGMA: val = bsq / rho; break;
              case BETA: val = 2.0 * p / bsq; break;
              default: val = 0.0;
              }
            }
            cellVal[cid(I, J, K)] = val;
          }
    }

    // statistics + transform
    double minPos = std::numeric_limits<double>::infinity();
    rawMin = std::numeric_limits<double>::infinity();
    rawMax = -rawMin;
    for (int k = 0; k < N[2]; k++)
      for (int j = 0; j < N[1]; j++)
        for (int i = i0; i < i1; i++) {
          const double v = cellVal[cid(i, j, k)];
          if (!std::isfinite(v)) { numNonFinite++; continue; }
          rawMin = std::min(rawMin, v);
          rawMax = std::max(rawMax, v);
          if (v > 0.0) minPos = std::min(minPos, v);
          else numNonPositive++;
        }
    if (logScale) {
      if (!std::isfinite(minPos)) throw std::runtime_error("no positive values for --log");
      logFloor = minPos;
    }
    for (int k = 0; k < N[2]; k++)
      for (int j = 0; j < N[1]; j++)
        for (int i = i0; i < i1; i++) {
          double &v = cellVal[cid(i, j, k)];
          if (!std::isfinite(v)) v = logScale ? logFloor : 0.0;
          if (logScale) v = std::log10(std::max(v, logFloor));
        }
    printf("variable %s (%s): raw range [%g, %g], %zu non-positive, %zu non-finite cells%s\n",
           opt.var.c_str(), logScale ? "log10" : "linear", rawMin, rawMax, numNonPositive, numNonFinite,
           logScale ? " (set to the smallest positive value)" : "");
  }

  void Converter::computeVertices()
  {
    vertPos.assign(numVerts, vec3f(0.f));
    for (int i = 0; i <= N[0]; i++) {
      const double r = re[i];
      for (int j = jlo; j <= jhi; j++) {
        const double st = std::sin(the[j]), ct = std::cos(the[j]);
        for (int k = 0; k < nKv; k++) {
          const double ph = phe[k];
          vertPos[vid(i, j, k)] = vec3f((float)(r * st * std::cos(ph)), (float)(r * st * std::sin(ph)), (float)(r * ct));
        }
      }
      if (northPole) vertPos[vid(i, 0, 0)] = vec3f(0.f, 0.f, (float)r);
      if (southPole) vertPos[vid(i, N[1], 0)] = vec3f(0.f, 0.f, (float)-r);
    }

    // vertex value = average over the kept cells sharing the vertex (each cell
    // counts once, also where its corners collapse on the axis)
    std::vector<double> sum(numVerts, 0.0);
    std::vector<uint32_t> cnt(numVerts, 0);
    for (int k = 0; k < N[2]; k++)
      for (int j = 0; j < N[1]; j++)
        for (int i = i0; i < i1; i++) {
          const double v = cellVal[cid(i, j, k)];
          size_t ids[8];
          int n = 0;
          for (int dk = 0; dk < 2; dk++)
            for (int dj = 0; dj < 2; dj++)
              for (int di = 0; di < 2; di++) {
                const size_t id = vid(i + di, j + dj, k + dk);
                bool dup = false;
                for (int q = 0; q < n; q++) dup = dup || ids[q] == id;
                if (!dup) ids[n++] = id;
              }
          for (int q = 0; q < n; q++) { sum[ids[q]] += v; cnt[ids[q]]++; }
        }
    vertVal.assign(numVerts, 0.f);
    for (size_t v = 0; v < numVerts; v++)
      vertVal[v] = cnt[v] ? (float)(sum[v] / cnt[v]) : 0.f;
  }

  // recursive bisection of the logical cell box along the direction with the
  // largest angular extent (x1 = ln r is the radial "angle" dr/r), balanced by
  // cell count; every part is a box of cells = a spherical shell sector
  void Converter::partition(const Box &b, int n, std::vector<Box> &out) const
  {
    if (n == 1) { out.push_back(b); return; }
    double ext[3];
    ext[0] = x1e[b.hi[0]] - x1e[b.lo[0]];
    ext[1] = the[b.hi[1]] - the[b.lo[1]];
    double maxSin = std::max(std::sin(the[b.lo[1]]), std::sin(the[b.hi[1]]));
    if (the[b.lo[1]] <= 0.5 * M_PI && the[b.hi[1]] >= 0.5 * M_PI) maxSin = 1.0;
    ext[2] = (phe[b.hi[2]] - phe[b.lo[2]]) * maxSin;
    int d = -1;
    for (int q = 0; q < 3; q++)
      if (b.hi[q] - b.lo[q] >= 2 && (d < 0 || ext[q] > ext[d])) d = q;
    if (d < 0) throw std::runtime_error("too many parts for the grid");
    const int nl = n / 2, nr = n - nl;
    const int len = b.hi[d] - b.lo[d];
    int cut = b.lo[d] + (int)std::lround((double)len * nl / n);
    cut = std::max(b.lo[d] + 1, std::min(b.hi[d] - 1, cut));
    Box l = b, r = b;
    l.hi[d] = cut;
    r.lo[d] = cut;
    partition(l, nl, out);
    partition(r, nr, out);
  }

  // Field at a Cartesian point, trilinear in the logical (x1, x2, x3) cell
  // coordinates of the per-vertex values (NaN outside the kept grid). Used by
  // the CPU reference renderer only.
  double Converter::sampleLogical(double x, double y, double z) const
  {
    const double r = std::sqrt(x * x + y * y + z * z);
    if (r - opt.mksR0 <= 0.0) return NAN;
    const double x1 = std::log(r - opt.mksR0);
    const double th = std::acos(std::max(-1.0, std::min(1.0, z / r)));
    double x2 = th;
    for (int it = 0; it < 8; it++)   // invert theta = x2 + h/2 sin(2 x2)
      x2 -= (x2 + 0.5 * opt.mksH * std::sin(2.0 * x2) - th) / (1.0 + opt.mksH * std::cos(2.0 * x2));
    double x3 = std::atan2(y, x);
    if (x3 < phe[0]) x3 += 2.0 * M_PI;
    const double u[3] = { (x1 - x1e[0]) / (x1e[N[0]] - x1e[0]) * N[0],
                          (x2 - x2e[0]) / (x2e[N[1]] - x2e[0]) * N[1],
                          (x3 - x3e[0]) / (x3e[N[2]] - x3e[0]) * N[2] };
    int c[3];
    double f[3];
    for (int d = 0; d < 3; d++) {
      if (u[d] < 0.0 || u[d] > N[d]) return NAN;
      c[d] = std::min(N[d] - 1, (int)u[d]);
      f[d] = u[d] - c[d];
    }
    if (c[0] < i0 || c[0] >= i1) return NAN;
    double v = 0.0;
    for (int dk = 0; dk < 2; dk++)
      for (int dj = 0; dj < 2; dj++)
        for (int di = 0; di < 2; di++)
          v += (di ? f[0] : 1 - f[0]) * (dj ? f[1] : 1 - f[1]) * (dk ? f[2] : 1 - f[2]) *
               vertVal[vid(c[0] + di, c[1] + dj, c[2] + dk)];
    return v;
  }

  // CPU reference: first crossing of the iso level along each primary ray,
  // refined by bisection, Lambert shading from the gradient; writes a PPM
  void Converter::refRender() const
  {
    const int W = opt.refRes[0], H = opt.refRes[1];
    const double *c = opt.refCam;
    double eye[3] = { c[0], c[1], c[2] }, dir[3] = { c[3] - c[0], c[4] - c[1], c[5] - c[2] };
    auto norm = [](double *v) { double l = std::sqrt(dot3(v, v)); for (int q = 0; q < 3; q++) v[q] /= l; };
    norm(dir);
    double up[3] = { c[6], c[7], c[8] }, right[3], camUp[3];
    cross3(dir, up, right); norm(right);
    cross3(right, dir, camUp); norm(camUp);
    const double th = std::tan(0.5 * opt.refFovy * M_PI / 180.0), aspect = (double)W / H;
    // light: from above-left of the camera
    double light[3];
    for (int q = 0; q < 3; q++) light[q] = -dir[q] + 0.6 * camUp[q] - 0.5 * right[q];
    norm(light);
    const double iso = opt.refIso, R = opt.refRmax, dt = opt.refStep;
    std::vector<unsigned char> img((size_t)W * H * 3, 0);
    size_t hits = 0;
#pragma omp parallel for schedule(dynamic, 4) reduction(+:hits)
    for (int py = 0; py < H; py++)
      for (int px = 0; px < W; px++) {
        const double sx = (2.0 * (px + 0.5) / W - 1.0) * th * aspect, sy = (1.0 - 2.0 * (py + 0.5) / H) * th;
        double d[3];
        for (int q = 0; q < 3; q++) d[q] = dir[q] + sx * right[q] + sy * camUp[q];
        norm(d);
        // ray / sphere r < R
        const double b = dot3(eye, d), cc = dot3(eye, eye) - R * R, disc = b * b - cc;
        if (disc <= 0) continue;
        const double t0 = std::max(0.0, -b - std::sqrt(disc)), t1 = -b + std::sqrt(disc);
        double tPrev = t0, fPrev = NAN;
        for (double t = t0; t <= t1; t += dt) {
          const double f = sampleLogical(eye[0] + t * d[0], eye[1] + t * d[1], eye[2] + t * d[2]);
          if (!std::isnan(f) && !std::isnan(fPrev) && fPrev < iso && f >= iso) {
            double a = tPrev, bb = t;
            for (int it = 0; it < 20; it++) {
              const double m = 0.5 * (a + bb);
              const double fm = sampleLogical(eye[0] + m * d[0], eye[1] + m * d[1], eye[2] + m * d[2]);
              if (!std::isnan(fm) && fm >= iso) bb = m; else a = m;
            }
            const double P[3] = { eye[0] + bb * d[0], eye[1] + bb * d[1], eye[2] + bb * d[2] };
            const double e = 0.02;
            double g[3];
            for (int q = 0; q < 3; q++) {
              double Pp[3] = { P[0], P[1], P[2] }, Pm[3] = { P[0], P[1], P[2] };
              Pp[q] += e; Pm[q] -= e;
              const double fp = sampleLogical(Pp[0], Pp[1], Pp[2]), fmn = sampleLogical(Pm[0], Pm[1], Pm[2]);
              g[q] = (std::isnan(fp) || std::isnan(fmn)) ? 0.0 : (fp - fmn);
            }
            double n[3] = { -g[0], -g[1], -g[2] };   // outward = towards lower values
            const double gl = std::sqrt(dot3(n, n));
            double shade = 0.5;
            if (gl > 0) { for (int q = 0; q < 3; q++) n[q] /= gl; shade = 0.15 + 0.85 * std::max(0.0, dot3(n, light)); }
            const double rP = std::sqrt(dot3(P, P));
            const double hue = std::min(1.0, std::max(0.0, (rP - 5.0) / 60.0));   // colour by radius
            unsigned char *pix = &img[((size_t)py * W + px) * 3];
            pix[0] = (unsigned char)std::min(255.0, 255.0 * shade * (0.95));
            pix[1] = (unsigned char)std::min(255.0, 255.0 * shade * (0.55 + 0.35 * hue));
            pix[2] = (unsigned char)std::min(255.0, 255.0 * shade * (0.35 + 0.5 * (1 - hue)));
            hits++;
            break;
          }
          tPrev = t;
          fPrev = f;
        }
      }
    FILE *fp = fopen(opt.refRender.c_str(), "wb");
    if (!fp) throw std::runtime_error("cannot write " + opt.refRender);
    fprintf(fp, "P6\n%d %d\n255\n", W, H);
    fwrite(img.data(), 1, img.size(), fp);
    fclose(fp);
    printf("reference render %s: %dx%d, iso %g, %zu pixels hit\n", opt.refRender.c_str(), W, H, iso, hits);
  }

  struct PartStats {
    std::string file;
    Box box;
    size_t hexes = 0, wedges = 0, verts = 0;
    umesh::box3f bounds;
    umesh::range1f range;
    size_t badHex = 0, badWedge = 0, tets = 0, badTet = 0;
    double volume = 0.0, minHexVol = INFINITY, minWedgeVol = INFINITY, minTetVol = INFINITY;
  };

  std::string jsonEscape(const std::string &s)
  {
    std::string o;
    for (char c : s) { if (c == '"' || c == '\\') o += '\\'; o += c; }
    return o;
  }

  void Converter::run()
  {
    setup();
    if (opt.listVars) return;
    computeCellValues();
    computeVertices();
    if (!opt.refRender.empty()) { refRender(); return; }

    Box all;
    all.lo[0] = i0; all.hi[0] = i1;
    all.lo[1] = 0; all.hi[1] = N[1];
    all.lo[2] = 0; all.hi[2] = N[2];
    std::vector<Box> boxes;
    partition(all, opt.numParts, boxes);

    std::string baseName = opt.outBase;
    const size_t slash = baseName.find_last_of('/');
    if (slash != std::string::npos) baseName = baseName.substr(slash + 1);

    std::vector<PartStats> stats(boxes.size());
    std::vector<int> g2l(numVerts, -1);
    std::vector<size_t> touched;
    for (size_t p = 0; p < boxes.size(); p++) {
      const Box &b = boxes[p];
      PartStats &st = stats[p];
      st.box = b;
      auto mesh = std::make_shared<UMesh>();
      mesh->perVertex = std::make_shared<umesh::Attribute>();
      mesh->perVertex->name = opt.var;
      touched.clear();
      auto L = [&](size_t gid) -> int {
        int &l = g2l[gid];
        if (l < 0) {
          l = (int)mesh->vertices.size();
          mesh->vertices.push_back(vertPos[gid]);
          mesh->perVertex->values.push_back(vertVal[gid]);
          touched.push_back(gid);
        }
        return l;
      };
      if (opt.tets) mesh->tets.reserve(6 * b.numCells()); else mesh->hexes.reserve(b.numCells());
      for (int k = b.lo[2]; k < b.hi[2]; k++)
        for (int j = b.lo[1]; j < b.hi[1]; j++)
          for (int i = b.lo[0]; i < b.hi[0]; i++) {
            if (opt.tets) {
              // Tet mode: the logical hex (VTK order, as below) split into the 6
              // tets around its diagonal 0-6. The pattern is the same in every
              // cell, so the face diagonals of neighbours match (no cracks), also
              // across the phi seam; on the axis the collapsed corners share one
              // vertex id and the 3 degenerate tets are dropped (-> 3 tets).
              const size_t g[8] = { vid(i, j, k), vid(i + 1, j, k), vid(i + 1, j + 1, k), vid(i, j + 1, k),
                                    vid(i, j, k + 1), vid(i + 1, j, k + 1), vid(i + 1, j + 1, k + 1),
                                    vid(i, j + 1, k + 1) };
              static const int T[6][4] = { {0,1,2,6}, {0,2,3,6}, {0,3,7,6}, {0,7,4,6}, {0,4,5,6}, {0,5,1,6} };
              for (auto &t : T) {
                if (g[t[0]] == g[t[1]] || g[t[0]] == g[t[2]] || g[t[0]] == g[t[3]] ||
                    g[t[1]] == g[t[2]] || g[t[1]] == g[t[3]] || g[t[2]] == g[t[3]])
                  continue;
                UMesh::Tet tet(L(g[t[0]]), L(g[t[1]]), L(g[t[2]]), L(g[t[3]]));
                // barney's tet test needs dot(v3-v0, cross(v1-v0, v2-v0)) > 0
                const double vol = tetVol6(mesh->vertices[tet.x], mesh->vertices[tet.y],
                                           mesh->vertices[tet.z], mesh->vertices[tet.w]) / 6.0;
                if (vol <= 0.0) st.badTet++;
                st.volume += vol; st.minTetVol = std::min(st.minTetVol, vol);
                mesh->tets.push_back(tet);
              }
              continue;
            }
            if (northPole && j == 0) {
              // Wedge in umesh/VTK order: front triangle (0,1,2) on the inner
              // sphere r_i, back triangle (3,4,5) on r_{i+1}, v3 above v0 etc.;
              // the right-hand normal of (0,1,2) points away from (3,4,5)
              // (-r direction), as required by fixNegativeVolumeElements
              int v[6] = { L(vid(i, 0, k)),     L(vid(i, 1, k + 1)),     L(vid(i, 1, k)),
                           L(vid(i + 1, 0, k)), L(vid(i + 1, 1, k + 1)), L(vid(i + 1, 1, k)) };
              UMesh::Wedge w(v[0], v[1], v[2], v[3], v[4], v[5]);
              vec3f P[6];
              for (int q = 0; q < 6; q++) P[q] = mesh->vertices[v[q]];
              if (!wedgeOrientedOK(P)) st.badWedge++;
              const double vol = wedgeVolume(P);
              st.volume += vol; st.minWedgeVol = std::min(st.minWedgeVol, vol);
              mesh->wedges.push_back(w);
            } else if (southPole && j == N[1] - 1) {
              int v[6] = { L(vid(i, N[1], k)),     L(vid(i, N[1] - 1, k)),     L(vid(i, N[1] - 1, k + 1)),
                           L(vid(i + 1, N[1], k)), L(vid(i + 1, N[1] - 1, k)), L(vid(i + 1, N[1] - 1, k + 1)) };
              UMesh::Wedge w(v[0], v[1], v[2], v[3], v[4], v[5]);
              vec3f P[6];
              for (int q = 0; q < 6; q++) P[q] = mesh->vertices[v[q]];
              if (!wedgeOrientedOK(P)) st.badWedge++;
              const double vol = wedgeVolume(P);
              st.volume += vol; st.minWedgeVol = std::min(st.minWedgeVol, vol);
              mesh->wedges.push_back(w);
            } else {
              // Hex in umesh/VTK order: base quad (0,1,2,3) = (i,j),(i+1,j),(i+1,j+1),(i,j+1)
              // at phi_k, top quad (4..7) the same at phi_{k+1}. (r, theta, phi) is
              // right handed, so the base normal points to the top (positive Jacobian)
              int v[8] = { L(vid(i, j, k)),         L(vid(i + 1, j, k)),
                           L(vid(i + 1, j + 1, k)), L(vid(i, j + 1, k)),
                           L(vid(i, j, k + 1)),     L(vid(i + 1, j, k + 1)),
                           L(vid(i + 1, j + 1, k + 1)), L(vid(i, j + 1, k + 1)) };
              UMesh::Hex h(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
              vec3f P[8];
              for (int q = 0; q < 8; q++) P[q] = mesh->vertices[v[q]];
              if (!hexOrientedOK(P)) st.badHex++;
              const double vol = hexVolume(P);
              st.volume += vol; st.minHexVol = std::min(st.minHexVol, vol);
              mesh->hexes.push_back(h);
            }
          }
      for (size_t g : touched) g2l[g] = -1;

      mesh->finalize();
      st.file = baseName + "_" + std::to_string(p) + ".umesh";
      mesh->saveTo(opt.outBase + "_" + std::to_string(p) + ".umesh");
      st.hexes = mesh->hexes.size();
      st.wedges = mesh->wedges.size();
      st.tets = mesh->tets.size();
      st.verts = mesh->vertices.size();
      st.bounds = mesh->getBounds();
      st.range = mesh->perVertex->valueRange;
      printf("part %3zu: cells i[%d,%d) j[%d,%d) k[%d,%d): %zu hexes, %zu wedges, %zu tets, %zu vertices, value [%g, %g]%s\n",
             p, b.lo[0], b.hi[0], b.lo[1], b.hi[1], b.lo[2], b.hi[2], st.hexes, st.wedges, st.tets, st.verts,
             st.range.lower, st.range.upper,
             (st.badHex || st.badWedge || st.badTet) ? "  ORIENTATION ERRORS" : "");
    }

    // ---- totals ----
    size_t hexes = 0, wedges = 0, verts = 0, badHex = 0, badWedge = 0, tets = 0, badTet = 0;
    double volume = 0.0, minHexVol = INFINITY, minWedgeVol = INFINITY, minTetVol = INFINITY;
    umesh::box3f bounds;
    umesh::range1f range;
    for (auto &s : stats) {
      hexes += s.hexes; wedges += s.wedges; verts += s.verts;
      badHex += s.badHex; badWedge += s.badWedge;
      tets += s.tets; badTet += s.badTet;
      minTetVol = std::min(minTetVol, s.minTetVol);
      volume += s.volume;
      minHexVol = std::min(minHexVol, s.minHexVol);
      minWedgeVol = std::min(minWedgeVol, s.minWedgeVol);
      bounds.extend(s.bounds);
      range.extend(s.range.lower); range.extend(s.range.upper);
    }
    const double shellVol = 4.0 / 3.0 * M_PI * (std::pow(re[i1], 3) - std::pow(re[i0], 3));

    // percentiles of the vertex values of the kept vertices (for choosing a
    // transfer-function range)
    std::vector<float> vv;
    vv.reserve(numVerts);
    for (int i = i0; i <= i1; i++)
      for (int j = 0; j <= N[1]; j++)
        for (int k = 0; k < nKv; k++) {
          if ((northPole && j == 0) || (southPole && j == N[1])) { if (k) continue; }
          vv.push_back(vertVal[vid(i, j, k)]);
        }
    std::sort(vv.begin(), vv.end());
    auto pct = [&](double q) { return vv.empty() ? 0.f : vv[std::min(vv.size() - 1, (size_t)(q * (vv.size() - 1) + 0.5))]; };

    printf("total: %zu parts, %zu elements = %zu hexes + %zu wedges + %zu tets, %zu vertices (%zu unique)\n",
           stats.size(), hexes + wedges + tets, hexes, wedges, tets, verts, vv.size());
    printf("bounds (%g %g %g) - (%g %g %g), value range [%g, %g]\n",
           bounds.lower.x, bounds.lower.y, bounds.lower.z, bounds.upper.x, bounds.upper.y, bounds.upper.z,
           range.lower, range.upper);
    printf("orientation (umesh fixNegativeVolumeElements criterion): %zu bad hexes, %zu bad wedges, %zu non-positive tets\n",
           badHex, badWedge, badTet);
    printf("signed volume: sum %.6g (spherical shell %.6g, ratio %.6f), min hex %.3g, min wedge %.3g, min tet %.3g\n",
           volume, shellVol, volume / shellVol, minHexVol, minWedgeVol, minTetVol);
    printf("vertex value percentiles: 1%% %g, 5%% %g, 50%% %g, 95%% %g, 99%% %g, 99.9%% %g\n",
           pct(0.01), pct(0.05), pct(0.5), pct(0.95), pct(0.99), pct(0.999));

    // ---- JSON summary ----
    std::ofstream js(opt.outBase + ".json");
    js.precision(9);
    js << "{\n";
    js << "  \"source\": \"" << jsonEscape(opt.datFile) << "\",\n";
    js << "  \"par\": \"" << jsonEscape(opt.parFile) << "\",\n";
    js << "  \"time\": " << snap.time << ",\n  \"iteration\": " << snap.iteration << ",\n";
    js << "  \"spin\": " << spin << ",\n  \"mass\": " << mass << ",\n  \"gamma\": " << gammaAd << ",\n";
    js << "  \"mks_h\": " << opt.mksH << ",\n  \"mks_R0\": " << opt.mksR0 << ",\n";
    js << "  \"grid_cells\": [" << N[0] << ", " << N[1] << ", " << N[2] << "],\n";
    js << "  \"kept_i_range\": [" << i0 << ", " << i1 << "],\n";
    js << "  \"r_range\": [" << re[i0] << ", " << re[i1] << "],\n";
    js << "  \"variable\": \"" << jsonEscape(opt.var) << "\",\n";
    js << "  \"scale\": \"" << (logScale ? "log10" : "linear") << "\",\n";
    js << "  \"raw_cell_range\": [" << rawMin << ", " << rawMax << "],\n";
    js << "  \"non_positive_cells\": " << numNonPositive << ",\n  \"non_finite_cells\": " << numNonFinite << ",\n";
    if (logScale) js << "  \"log_floor\": " << logFloor << ",\n";
    js << "  \"vertex_order\": \"umesh/VTK: hex base(0-3) at phi_k -> top(4-7) at phi_k+1, positive Jacobian; "
          "wedge front(0-2) on r_i, back(3-5) on r_i+1, front normal pointing away from back\",\n";
    js << "  \"num_parts\": " << stats.size() << ",\n";
    js << "  \"elements\": " << hexes + wedges + tets << ",\n  \"hexes\": " << hexes << ",\n  \"wedges\": " << wedges
       << ",\n  \"tets\": " << tets << ",\n";
    js << "  \"vertices\": " << verts << ",\n  \"unique_vertices\": " << vv.size() << ",\n";
    js << "  \"bounds\": [[" << bounds.lower.x << ", " << bounds.lower.y << ", " << bounds.lower.z << "], ["
       << bounds.upper.x << ", " << bounds.upper.y << ", " << bounds.upper.z << "]],\n";
    js << "  \"value_range\": [" << range.lower << ", " << range.upper << "],\n";
    js << "  \"value_percentiles\": {\"1\": " << pct(0.01) << ", \"5\": " << pct(0.05) << ", \"50\": " << pct(0.5)
       << ", \"95\": " << pct(0.95) << ", \"99\": " << pct(0.99) << ", \"99.9\": " << pct(0.999) << "},\n";
    js << "  \"bad_orientation\": {\"hexes\": " << badHex << ", \"wedges\": " << badWedge << "},\n";
    js << "  \"volume\": " << volume << ",\n  \"shell_volume\": " << shellVol << ",\n";
    js << "  \"parts\": [\n";
    for (size_t p = 0; p < stats.size(); p++) {
      const PartStats &s = stats[p];
      js << "    {\"file\": \"" << jsonEscape(s.file) << "\", \"ijk_lo\": [" << s.box.lo[0] << ", " << s.box.lo[1] << ", "
         << s.box.lo[2] << "], \"ijk_hi\": [" << s.box.hi[0] << ", " << s.box.hi[1] << ", " << s.box.hi[2] << "], "
         << "\"elements\": " << s.hexes + s.wedges + s.tets << ", \"hexes\": " << s.hexes << ", \"wedges\": " << s.wedges
         << ", \"tets\": " << s.tets
         << ", \"vertices\": " << s.verts << ", \"bounds\": [[" << s.bounds.lower.x << ", " << s.bounds.lower.y << ", "
         << s.bounds.lower.z << "], [" << s.bounds.upper.x << ", " << s.bounds.upper.y << ", " << s.bounds.upper.z
         << "]], \"value_range\": [" << s.range.lower << ", " << s.range.upper << "]}"
         << (p + 1 < stats.size() ? "," : "") << "\n";
    }
    js << "  ]\n}\n";
    if (!js) throw std::runtime_error("cannot write " + opt.outBase + ".json");

    if (badHex || badWedge || badTet)
      throw std::runtime_error("elements with wrong orientation were written");
  }

} // namespace

int main(int ac, char **av)
{
  try {
    Converter c;
    c.opt = parseArgs(ac, av);
    c.run();
  } catch (const std::exception &e) {
    std::cerr << "bhac2umesh: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
