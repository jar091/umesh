// foam2umesh - per-processor conversion (see Convert.h for conventions)

#include "Convert.h"
#include "FoamIO.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <set>
#include <cstring>
#include <cstdio>

namespace f2u {

  void ProcStats::add(const ProcStats &o)
  {
    cells += o.cells; points += o.points; faces += o.faces;
    internalFaces += o.internalFaces;
    for (int i = 0; i < 3; i++) facesBySize[i] += o.facesBySize[i];
    nativeHex += o.nativeHex; nativeWedge += o.nativeWedge;
    nativePyr += o.nativePyr; nativeTet += o.nativeTet;
    decomposedPoly += o.decomposedPoly;
    decomposedQuality += o.decomposedQuality;
    decomposedTopology += o.decomposedTopology;
    outVertices += o.outVertices;
    addedCellCentres += o.addedCellCentres;
    addedFaceCentres += o.addedFaceCentres;
    outTets += o.outTets; outPyrs += o.outPyrs;
    outWedges += o.outWedges; outHexes += o.outHexes;
    invertedTets += o.invertedTets; invertedPyrs += o.invertedPyrs;
    invertedWedges += o.invertedWedges; invertedHexes += o.invertedHexes;
    notStarTets += o.notStarTets; notStarPyrs += o.notStarPyrs;
    notStarWedges += o.notStarWedges; notStarHexes += o.notStarHexes;
    badCornerHexes += o.badCornerHexes; badCornerWedges += o.badCornerWedges;
    quadsSplitToTets += o.quadsSplitToTets;
    seamPoints += o.seamPoints; seamMatched += o.seamMatched;
    seamContributions += o.seamContributions;
    maskedCells += o.maskedCells; maskedVertices += o.maskedVertices;
    nanCells += o.nanCells;
    cellMin = std::min(cellMin, o.cellMin);
    cellMax = std::max(cellMax, o.cellMax);
  }

  size_t ProcMesh::memoryBytes() const
  {
    return vertices.capacity() * sizeof(vertices[0])
      + values.capacity() * sizeof(float)
      + tets.capacity() * sizeof(umesh::Tet)
      + pyrs.capacity() * sizeof(umesh::Pyr)
      + wedges.capacity() * sizeof(umesh::Wedge)
      + hexes.capacity() * sizeof(umesh::Hex);
  }

  namespace {

    using V3 = std::array<double,3>;

    inline V3 sub(const V3 &a, const V3 &b) { return { a[0]-b[0], a[1]-b[1], a[2]-b[2] }; }
    inline V3 add(const V3 &a, const V3 &b) { return { a[0]+b[0], a[1]+b[1], a[2]+b[2] }; }
    inline V3 mul(double s, const V3 &a) { return { s*a[0], s*a[1], s*a[2] }; }
    inline V3 cross(const V3 &a, const V3 &b)
    { return { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] }; }
    inline double dot(const V3 &a, const V3 &b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
    inline double len(const V3 &a) { return std::sqrt(dot(a,a)); }

    /*! all state needed while converting one processor */
    struct Converter {
      const ConvertOptions &opt;
      ProcMesh             &out;
      ProcStats            &st;

      std::vector<double>  P;          // points, xyz
      std::vector<int32_t> fOfs, fVtx; // faces (compact)
      std::vector<int32_t> owner, nei;
      std::vector<double>  cellVal;
      int64_t nP = 0, nF = 0, nIF = 0, nC = 0;

      std::vector<int32_t> cStart, cFaces; // cell -> faces (2*f + isNeighbour)

      std::vector<double>  vSum;           // per original point
      std::vector<int32_t> vCnt;

      std::vector<int32_t> faceCentreIdx;  // per face, -1 if none yet
      std::vector<std::pair<int32_t,int32_t>> pendingFaceCentres; // (vertex, face)

      // per-cell scratch
      std::vector<int32_t> fb;     // outward face vertex lists, concatenated
      std::vector<int32_t> fbOfs;  // offsets into fb (numFaces+1)
      std::vector<int32_t> fbFace; // original face id per local face
      std::vector<int32_t> uv;     // unique vertices of the cell

      Converter(const ConvertOptions &o, ProcMesh &m) : opt(o), out(m), st(m.stats) {}

      V3 pt(int32_t i) const { return { P[3*i], P[3*i+1], P[3*i+2] }; }

      V3 vtx(int32_t i) const {
        if (i < nP) return pt(i);
        const umesh::vec3f &v = out.vertices[i];
        return { v.x, v.y, v.z };
      }

      int32_t addVertex(const V3 &p, float value)
      {
        const int32_t idx = (int32_t)out.vertices.size();
        out.vertices.push_back(umesh::vec3f((float)p[0], (float)p[1], (float)p[2]));
        out.values.push_back(value);
        return idx;
      }

      int nFacesOfCell() const { return (int)fbOfs.size() - 1; }
      int faceSize(int lf) const { return fbOfs[lf+1] - fbOfs[lf]; }
      const int32_t *face(int lf) const { return fb.data() + fbOfs[lf]; }

      /*! gather outward-oriented faces of cell c into fb/fbOfs */
      void gatherCell(int64_t c)
      {
        fb.clear(); fbOfs.clear(); fbFace.clear(); uv.clear();
        fbOfs.push_back(0);
        for (int32_t k = cStart[c]; k < cStart[c+1]; k++) {
          const int32_t code = cFaces[k];
          const int32_t f = code >> 1;
          const bool flip = (code & 1) != 0;
          const int32_t b = fOfs[f], n = fOfs[f+1] - fOfs[f];
          if (!flip) {
            for (int i = 0; i < n; i++) fb.push_back(fVtx[b+i]);
          } else {
            // reversed orientation (first vertex kept), normal now
            // points out of the neighbour cell
            fb.push_back(fVtx[b]);
            for (int i = n-1; i >= 1; i--) fb.push_back(fVtx[b+i]);
          }
          fbOfs.push_back((int32_t)fb.size());
          fbFace.push_back(f);
        }
        uv.assign(fb.begin(), fb.end());
        std::sort(uv.begin(), uv.end());
        uv.erase(std::unique(uv.begin(), uv.end()), uv.end());
      }

      static int indexIn(const int32_t *a, int n, int32_t v)
      {
        for (int i = 0; i < n; i++) if (a[i] == v) return i;
        return -1;
      }

      /*! for every vertex of 'base' (n verts), find the vertex connected
          to it by the lateral edge in the given side faces. Returns false
          on inconsistent topology. */
      bool findPartners(const int32_t *base, int n, int baseFace, int otherFace,
                        int32_t *partner)
      {
        for (int i = 0; i < n; i++) partner[i] = -1;
        for (int lf = 0; lf < nFacesOfCell(); lf++) {
          if (lf == baseFace || lf == otherFace) continue;
          const int32_t *w = face(lf);
          const int m = faceSize(lf);
          if (m != 4) return false;
          int numBase = 0;
          for (int k = 0; k < m; k++) if (indexIn(base, n, w[k]) >= 0) numBase++;
          if (numBase != 2) return false;
          for (int k = 0; k < m; k++) {
            const int idx = indexIn(base, n, w[k]);
            if (idx < 0) continue;
            const int32_t prev = w[(k+m-1)%m], next = w[(k+1)%m];
            const bool inPrev = indexIn(base, n, prev) >= 0;
            const bool inNext = indexIn(base, n, next) >= 0;
            if (inPrev == inNext) return false;
            const int32_t cand = inPrev ? next : prev;
            if (partner[idx] < 0) partner[idx] = cand;
            else if (partner[idx] != cand) return false;
          }
        }
        for (int i = 0; i < n; i++) {
          if (partner[i] < 0) return false;
          if (indexIn(base, n, partner[i]) >= 0) return false;
          for (int j = 0; j < i; j++) if (partner[j] == partner[i]) return false;
        }
        // partners must form the opposite face
        const int32_t *o = face(otherFace);
        if (faceSize(otherFace) != n) return false;
        for (int i = 0; i < n; i++) if (indexIn(o, n, partner[i]) < 0) return false;
        return true;
      }

      bool tryHex(umesh::Hex &hex)
      {
        if (uv.size() != 8) return false;
        const int32_t *o = face(0);
        const int32_t b[4] = { o[0], o[3], o[2], o[1] }; // inward normal
        // opposite face: the one sharing no vertex with face 0
        int other = -1;
        for (int lf = 1; lf < 6; lf++) {
          bool shares = false;
          for (int k = 0; k < 4; k++) if (indexIn(b, 4, face(lf)[k]) >= 0) shares = true;
          if (!shares) { if (other >= 0) return false; other = lf; }
        }
        if (other < 0) return false;
        int32_t t[4];
        if (!findPartners(b, 4, 0, other, t)) return false;
        hex = umesh::Hex(b[0],b[1],b[2],b[3],t[0],t[1],t[2],t[3]);
        return true;
      }

      bool tryWedge(umesh::Wedge &w)
      {
        if (uv.size() != 6) return false;
        int tri[2], nt = 0;
        for (int lf = 0; lf < 5; lf++) if (faceSize(lf) == 3) tri[nt++] = lf;
        if (nt != 2) return false;
        const int32_t *a = face(tri[0]);
        // (v0,v1,v2) = outward oriented triangle (VTK wedge convention)
        const int32_t b[3] = { a[0], a[1], a[2] };
        int32_t t[3];
        if (!findPartners(b, 3, tri[0], tri[1], t)) return false;
        w = umesh::Wedge(b[0],b[1],b[2],t[0],t[1],t[2]);
        return true;
      }

      bool tryPyr(umesh::Pyr &p)
      {
        if (uv.size() != 5) return false;
        int q = -1;
        for (int lf = 0; lf < 5; lf++) if (faceSize(lf) == 4) q = lf;
        if (q < 0) return false;
        const int32_t *o = face(q);
        const int32_t b[4] = { o[0], o[3], o[2], o[1] }; // inward normal -> apex
        int32_t apex = -1;
        for (int32_t v : uv) if (indexIn(b, 4, v) < 0) apex = v;
        if (apex < 0) return false;
        for (int lf = 0; lf < 5; lf++) {
          if (lf == q) continue;
          if (faceSize(lf) != 3 || indexIn(face(lf), 3, apex) < 0) return false;
        }
        p = umesh::Pyr(b[0],b[1],b[2],b[3],apex);
        return true;
      }

      bool tryTet(umesh::Tet &t)
      {
        if (uv.size() != 4) return false;
        const int32_t *a = face(0);
        int32_t d = -1;
        for (int32_t v : uv) if (indexIn(a, 3, v) < 0) d = v;
        if (d < 0) return false;
        t = umesh::Tet(a[0], a[2], a[1], d); // outward (a,b,c) -> (a,c,b)
        return true;
      }

      template<int N>
      void coords(const int32_t *idx, V3 *v) const
      { for (int i = 0; i < N; i++) v[i] = vtx(idx[i]); }

      /*! face centre (area weighted, as OpenFOAM) and area vector of an
          outward face given by its vertex list */
      void faceGeom(const int32_t *w, int n, V3 &fc, V3 &Sf) const
      {
        if (n == 3) {
          const V3 a = vtx(w[0]), b = vtx(w[1]), c = vtx(w[2]);
          fc = mul(1.0/3.0, add(add(a,b),c));
          Sf = mul(0.5, cross(sub(b,a), sub(c,a)));
          return;
        }
        V3 avg = {0,0,0};
        for (int i = 0; i < n; i++) avg = add(avg, vtx(w[i]));
        avg = mul(1.0/n, avg);
        V3 sumN = {0,0,0}, sumAc = {0,0,0};
        double sumA = 0.0;
        for (int i = 0; i < n; i++) {
          const V3 a = vtx(w[i]), b = vtx(w[(i+1)%n]);
          const V3 c = add(add(a,b),avg);
          const V3 nn = cross(sub(b,a), sub(avg,a));
          const double ar = len(nn);
          sumN = add(sumN, nn);
          sumA += ar;
          sumAc = add(sumAc, mul(ar, c));
        }
        fc = (sumA > 1e-300) ? mul(1.0/(3.0*sumA), sumAc) : avg;
        Sf = mul(0.5, sumN);
      }

      /*! decompose the current cell into pyramids (quad faces) and tets
          (triangle faces; faces with >4 vertices are fanned around a
          shared face-centre vertex) towards an added cell-centre vertex */
      void decompose(float value)
      {
        const int nf = nFacesOfCell();
        // --- cell centre as in OpenFOAM (primitiveMeshCellCentresAndVols)
        V3 cEst = {0,0,0};
        std::vector<V3> fcs(nf), sfs(nf);
        for (int lf = 0; lf < nf; lf++) {
          faceGeom(face(lf), faceSize(lf), fcs[lf], sfs[lf]);
          cEst = add(cEst, fcs[lf]);
        }
        cEst = mul(1.0/nf, cEst);
        V3 sumVc = {0,0,0};
        double sumV = 0.0;
        for (int lf = 0; lf < nf; lf++) {
          const double pv = std::max(dot(sfs[lf], sub(fcs[lf], cEst)), 1e-300);
          const V3 pc = add(mul(0.75, fcs[lf]), mul(0.25, cEst));
          sumV += pv;
          sumVc = add(sumVc, mul(pv, pc));
        }
        const V3 C = (sumV > 1e-250) ? mul(1.0/sumV, sumVc) : cEst;
        const int32_t ci = addVertex(C, value);
        st.addedCellCentres++;

        for (int lf = 0; lf < nf; lf++) {
          const int32_t *w = face(lf);
          const int n = faceSize(lf);
          if (n == 3) {
            const int32_t e[4] = { w[0], w[2], w[1], ci };
            emitTet(e);
          } else if (n == 4) {
            const int32_t e[5] = { w[0], w[3], w[2], w[1], ci };
            emitDecompPyr(e);
          } else {
            const int32_t f = fbFace[lf];
            int32_t fci = faceCentreIdx[f];
            if (fci < 0) {
              fci = addVertex(fcs[lf], 0.f); // value filled in later
              faceCentreIdx[f] = fci;
              pendingFaceCentres.push_back({fci, f});
              st.addedFaceCentres++;
            }
            for (int i = 0; i < n; i++) {
              const int32_t e[4] = { fci, w[(i+1)%n], w[i], ci };
              emitTet(e);
            }
          }
        }
      }

      void emitTet(const int32_t *e)
      {
        V3 v[4]; coords<4>(e, v);
        if (!(tetVol(v[0],v[1],v[2],v[3]) > 0.0)) { st.invertedTets++; st.notStarTets++; }
        out.tets.push_back(umesh::Tet(e[0],e[1],e[2],e[3]));
        st.outTets++;
      }
      /*! pyramid of a decomposed cell: base (e0..e3), apex e4 = cell centre */
      void emitDecompPyr(const int32_t *e)
      {
        V3 v[5]; coords<5>(e, v);
        const VolCheck vc = checkElement(v, pyrFaces);
        if (!vc.allPositive && opt.splitBadPyramids) {
          // try both diagonals of the (warped) base quad
          for (int d = 0; d < 2; d++) {
            const int a = d, b = d+1, c = d+2, dd = (d+3)%4;
            const double v1 = tetVol(v[a],v[b],v[c],v[4]);
            const double v2 = tetVol(v[a],v[c],v[dd],v[4]);
            if (v1 > 0.0 && v2 > 0.0) {
              out.tets.push_back(umesh::Tet(e[a],e[b],e[c],e[4]));
              out.tets.push_back(umesh::Tet(e[a],e[c],e[dd],e[4]));
              st.outTets += 2;
              st.quadsSplitToTets++;
              return;
            }
          }
        }
        emitPyr(e, vc);
      }
      void emitPyr(const int32_t *e, const VolCheck &vc)
      {
        if (!(vc.total > 0.0)) st.invertedPyrs++;
        if (!vc.allPositive) st.notStarPyrs++;
        out.pyrs.push_back(umesh::Pyr(e[0],e[1],e[2],e[3],e[4]));
        st.outPyrs++;
      }
      void emitWedge(const umesh::Wedge &w, const VolCheck &vc)
      {
        if (!(vc.total > 0.0)) st.invertedWedges++;
        if (!vc.allPositive) st.notStarWedges++;
        out.wedges.push_back(w);
        st.outWedges++;
      }
      void emitHex(const umesh::Hex &h, const VolCheck &vc)
      {
        if (!(vc.total > 0.0)) st.invertedHexes++;
        if (!vc.allPositive) st.notStarHexes++;
        out.hexes.push_back(h);
        st.outHexes++;
      }

      bool keepNative(const VolCheck &vc, bool cornersOK) const
      {
        if (opt.qualityMode == 0) return true;
        if (!vc.allPositive) return false;
        if (opt.qualityMode == 2 && !cornersOK) return false;
        return true;
      }

      std::vector<uint8_t> cellMasked;  // per cell
      std::vector<uint8_t> vMask;       // per output vertex

      /*! cells within opt.wallMaskLayers layers of a wall patch */
      void computeWallMask()
      {
        cellMasked.assign(nC, 0);
        std::vector<int32_t> front, next;
        for (const auto &pa : patches) {
          if (pa.type != "wall") continue;
          for (int64_t f = pa.startFace; f < pa.startFace + pa.nFaces; f++) {
            const int32_t c = owner[f];
            if (!cellMasked[c]) { cellMasked[c] = 1; front.push_back(c); }
          }
        }
        for (int layer = 2; layer <= opt.wallMaskLayers && !front.empty(); layer++) {
          next.clear();
          for (int32_t c : front)
            for (int32_t k = cStart[c]; k < cStart[c+1]; k++) {
              const int32_t f = cFaces[k] >> 1;
              if (f >= nIF) continue; // boundary face: no neighbour here
              const int32_t o = (owner[f] == c) ? nei[f] : owner[f];
              if (!cellMasked[o]) { cellMasked[o] = 1; next.push_back(o); }
            }
          front.swap(next);
        }
        for (int64_t c = 0; c < nC; c++) st.maskedCells += cellMasked[c];
      }

      void convertCell(int64_t c)
      {
        const size_t t0 = out.tets.size(), p0 = out.pyrs.size(),
          w0 = out.wedges.size(), h0 = out.hexes.size();
        convertCellElements(c);
        if (!cellMasked.empty() && cellMasked[c]) {
          if (vMask.size() < out.vertices.size()) vMask.resize(out.vertices.size() + 1024, 0);
          for (size_t i = t0; i < out.tets.size(); i++) for (int k = 0; k < 4; k++) vMask[out.tets[i][k]] = 1;
          for (size_t i = p0; i < out.pyrs.size(); i++) for (int k = 0; k < 5; k++) vMask[out.pyrs[i][k]] = 1;
          for (size_t i = w0; i < out.wedges.size(); i++) for (int k = 0; k < 6; k++) vMask[out.wedges[i][k]] = 1;
          for (size_t i = h0; i < out.hexes.size(); i++) for (int k = 0; k < 8; k++) vMask[out.hexes[i][k]] = 1;
        }
      }

      void convertCellElements(int64_t c)
      {
        gatherCell(c);
        const int nf = nFacesOfCell();
        int n3 = 0, n4 = 0, nx = 0;
        for (int lf = 0; lf < nf; lf++) {
          const int s = faceSize(lf);
          if (s == 3) n3++; else if (s == 4) n4++; else nx++;
        }
        const float value = (float)cellVal[c];
        for (int32_t v : uv) { vSum[v] += cellVal[c]; vCnt[v]++; }

        bool shapeLike = false;
        if (nx == 0) {
          if (nf == 6 && n4 == 6) {
            shapeLike = true;
            umesh::Hex h;
            if (tryHex(h)) {
              V3 v[8]; coords<8>(&h.base.x, v);
              const VolCheck vc = checkElement(v, hexFaces);
              const bool corners = hexCornersOK(v);
              if (!corners) st.badCornerHexes++;
              if (keepNative(vc, corners)) {
                emitHex(h, vc); st.nativeHex++;
                return;
              }
              st.decomposedQuality++;
              decompose(value);
              return;
            }
          } else if (nf == 5 && n3 == 2 && n4 == 3) {
            shapeLike = true;
            umesh::Wedge w;
            if (tryWedge(w)) {
              V3 v[6]; coords<6>(&w.front.x, v);
              const VolCheck vc = checkElement(v, wedgeFaces);
              const bool corners = wedgeCornersOK(v);
              if (!corners) st.badCornerWedges++;
              if (keepNative(vc, corners)) {
                emitWedge(w, vc); st.nativeWedge++;
                return;
              }
              st.decomposedQuality++;
              decompose(value);
              return;
            }
          } else if (nf == 5 && n3 == 4 && n4 == 1) {
            shapeLike = true;
            umesh::Pyr p;
            if (tryPyr(p)) {
              V3 v[5]; coords<5>(&p.base.x, v);
              const VolCheck vc = checkElement(v, pyrFaces);
              if (keepNative(vc, true)) {
                emitPyr(&p.base.x, vc); st.nativePyr++;
                return;
              }
              st.decomposedQuality++;
              decompose(value);
              return;
            }
          } else if (nf == 4 && n3 == 4) {
            shapeLike = true;
            umesh::Tet t;
            if (tryTet(t)) {
              emitTet(&t.x); st.nativeTet++; // cannot be decomposed further
              return;
            }
          }
        }
        if (shapeLike) st.decomposedTopology++;
        else           st.decomposedPoly++;
        decompose(value);
      }

      std::vector<foam::Patch> patches;

      /*! read mesh, boundary and field of one processor, build cell->faces */
      void load(const std::string &caseDir, int procID)
      {
        const std::string pdir = caseDir + "/processor" + std::to_string(procID);
        const std::string mdir = pdir + "/constant/polyMesh/";
        out.procID = procID;

        foam::readPoints(mdir + "points", P);
        foam::readFaces(mdir + "faces", fOfs, fVtx);
        foam::Header ownerHdr;
        foam::readLabelList(mdir + "owner", owner, &ownerHdr);
        foam::readLabelList(mdir + "neighbour", nei);

        nP  = (int64_t)P.size() / 3;
        nF  = (int64_t)fOfs.size() - 1;
        nIF = (int64_t)nei.size();
        auto note = foam::parseNote(ownerHdr.note);
        if (note.count("nCells")) nC = note["nCells"];
        else {
          for (int32_t o : owner) nC = std::max<int64_t>(nC, o + 1);
          for (int32_t o : nei)   nC = std::max<int64_t>(nC, o + 1);
        }
        if ((int64_t)owner.size() != nF)
          throw std::runtime_error(mdir + ": owner size != number of faces");
        if (note.count("nPoints") && note["nPoints"] != nP)
          throw std::runtime_error(mdir + ": nPoints mismatch");
        if (note.count("nInternalFaces") && note["nInternalFaces"] != nIF)
          throw std::runtime_error(mdir + ": nInternalFaces mismatch");

        // --- scalar field
        {
          std::vector<double> raw;
          int nc = 0;
          const std::string fpath = pdir + "/" + opt.time + "/" + opt.field;
          foam::readInternalField(fpath, raw, nc);
          if ((int64_t)raw.size() != nc * nC)
            throw std::runtime_error(fpath + ": field size " + std::to_string(raw.size()/nc)
                                     + " != nCells " + std::to_string(nC));
          cellVal.resize(nC);
          for (int64_t c = 0; c < nC; c++) {
            double v;
            if (nc == 1) v = raw[c];
            else v = std::sqrt(raw[3*c]*raw[3*c] + raw[3*c+1]*raw[3*c+1] + raw[3*c+2]*raw[3*c+2]);
            if (!std::isfinite(v)) { st.nanCells++; v = 0.0; }
            cellVal[c] = v;
            st.cellMin = std::min(st.cellMin, v);
            st.cellMax = std::max(st.cellMax, v);
          }
        }

        st.cells = nC; st.points = nP; st.faces = nF; st.internalFaces = nIF;
        for (int64_t f = 0; f < nF; f++) {
          const int n = fOfs[f+1] - fOfs[f];
          st.facesBySize[n == 3 ? 0 : (n == 4 ? 1 : 2)]++;
        }

        // --- cell -> faces
        cStart.assign(nC + 1, 0);
        for (int64_t f = 0; f < nF; f++) cStart[owner[f] + 1]++;
        for (int64_t f = 0; f < nIF; f++) cStart[nei[f] + 1]++;
        for (int64_t c = 0; c < nC; c++) cStart[c+1] += cStart[c];
        cFaces.resize(cStart[nC]);
        {
          std::vector<int32_t> fill(cStart.begin(), cStart.end() - 1);
          for (int64_t f = 0; f < nF; f++) cFaces[fill[owner[f]]++] = int32_t(2*f);
          for (int64_t f = 0; f < nIF; f++) cFaces[fill[nei[f]]++]  = int32_t(2*f + 1);
        }

        foam::readBoundary(mdir + "boundary", patches);
      }

      /*! the unique points on processor patches, and the neighbouring
          processors */
      void procPatchPoints(std::vector<int32_t> &pts, std::vector<int32_t> &neighbours) const
      {
        pts.clear(); neighbours.clear();
        for (const auto &pa : patches) {
          if (pa.type != "processor" || pa.nFaces <= 0) continue;
          if (pa.neighbProcNo >= 0) neighbours.push_back(pa.neighbProcNo);
          for (int64_t f = pa.startFace; f < pa.startFace + pa.nFaces; f++)
            for (int k = fOfs[f]; k < fOfs[f+1]; k++) pts.push_back(fVtx[k]);
        }
        std::sort(pts.begin(), pts.end());
        pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
        std::sort(neighbours.begin(), neighbours.end());
        neighbours.erase(std::unique(neighbours.begin(), neighbours.end()), neighbours.end());
      }

      /*! accumulate incident cell values of every point (seam pre-pass) */
      void accumulateOnly()
      {
        vSum.assign(nP, 0.0);
        vCnt.assign(nP, 0);
        for (int64_t c = 0; c < nC; c++) {
          gatherCell(c);
          for (int32_t v : uv) { vSum[v] += cellVal[c]; vCnt[v]++; }
        }
      }

      struct SeamRecord { double x, y, z, sum; int64_t count; };
      static constexpr uint64_t seamMagic = 0x314d41455355324full; // "O2USEAM1"

      static std::string seamFile(const std::string &dir, int proc)
      { return dir + "/proc" + std::to_string(proc) + ".seam"; }

      void writeSeam(const std::string &dir)
      {
        std::vector<int32_t> pts, neighbours;
        procPatchPoints(pts, neighbours);
        const std::string fn = seamFile(dir, out.procID);
        const std::string tmp = fn + ".tmp";
        FILE *fp = fopen(tmp.c_str(), "wb");
        if (!fp) throw std::runtime_error("cannot write " + tmp);
        const int32_t id = out.procID, nn = (int32_t)neighbours.size();
        const int64_t np = (int64_t)pts.size();
        fwrite(&seamMagic, sizeof(seamMagic), 1, fp);
        fwrite(&id, sizeof(id), 1, fp);
        fwrite(&nn, sizeof(nn), 1, fp);
        fwrite(neighbours.data(), sizeof(int32_t), nn, fp);
        fwrite(&np, sizeof(np), 1, fp);
        std::vector<SeamRecord> recs(np);
        for (int64_t i = 0; i < np; i++) {
          const int32_t v = pts[i];
          recs[i] = { P[3*v], P[3*v+1], P[3*v+2], vSum[v], (int64_t)vCnt[v] };
        }
        fwrite(recs.data(), sizeof(SeamRecord), np, fp);
        fclose(fp);
        if (rename(tmp.c_str(), fn.c_str()) != 0)
          throw std::runtime_error("cannot rename " + tmp);
        st.seamPoints = np;
      }

      static void readSeam(const std::string &fn, std::vector<int32_t> &neighbours,
                           std::vector<SeamRecord> *recs)
      {
        FILE *fp = fopen(fn.c_str(), "rb");
        if (!fp) throw std::runtime_error("cannot open seam file " + fn
                                          + " (run the --seam-prepass first)");
        uint64_t magic = 0; int32_t id = 0, nn = 0; int64_t np = 0;
        bool ok = fread(&magic, sizeof(magic), 1, fp) == 1 && magic == seamMagic
          && fread(&id, sizeof(id), 1, fp) == 1 && fread(&nn, sizeof(nn), 1, fp) == 1;
        if (ok) {
          neighbours.resize(nn);
          ok = fread(neighbours.data(), sizeof(int32_t), nn, fp) == (size_t)nn
            && fread(&np, sizeof(np), 1, fp) == 1;
        }
        if (ok && recs) {
          recs->resize(np);
          ok = fread(recs->data(), sizeof(SeamRecord), np, fp) == (size_t)np;
        }
        fclose(fp);
        if (!ok) throw std::runtime_error("corrupt seam file " + fn);
      }

      struct BitKey {
        uint64_t a, b, c;
        bool operator==(const BitKey &o) const { return a == o.a && b == o.b && c == o.c; }
      };
      struct BitKeyHash {
        size_t operator()(const BitKey &k) const {
          uint64_t h = k.a * 0x9E3779B97F4A7C15ull;
          h ^= k.b + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
          h ^= k.c + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
          return (size_t)h;
        }
      };
      static BitKey key(double x, double y, double z)
      {
        BitKey k;
        memcpy(&k.a, &x, 8); memcpy(&k.b, &y, 8); memcpy(&k.c, &z, 8);
        return k;
      }

      /*! add the contributions of other processors sharing points on
          this processor's processor patches (matched by bitwise
          identical coordinates) */
      void applySeams(const std::string &dir)
      {
        std::vector<int32_t> pts, n1;
        procPatchPoints(pts, n1);
        st.seamPoints = (int64_t)pts.size();
        if (pts.empty()) return;
        // 2-ring of processors: face neighbours and their face neighbours
        std::set<int32_t> ring(n1.begin(), n1.end());
        for (int32_t q : n1) {
          std::vector<int32_t> n2;
          readSeam(seamFile(dir, q), n2, nullptr);
          ring.insert(n2.begin(), n2.end());
        }
        ring.erase(out.procID);

        std::unordered_map<BitKey, int32_t, BitKeyHash> map;
        map.reserve(pts.size() * 2);
        for (int32_t v : pts) map[key(P[3*v], P[3*v+1], P[3*v+2])] = v;
        std::vector<double>  addSum(pts.size(), 0.0);
        std::vector<int32_t> addCnt(pts.size(), 0);
        std::unordered_map<int32_t, int32_t> local; // point -> index in pts
        for (size_t i = 0; i < pts.size(); i++) local[pts[i]] = (int32_t)i;

        std::vector<SeamRecord> recs;
        std::vector<int32_t> dummy;
        for (int32_t r : ring) {
          readSeam(seamFile(dir, r), dummy, &recs);
          for (const auto &rec : recs) {
            auto it = map.find(key(rec.x, rec.y, rec.z));
            if (it == map.end()) continue;
            const int32_t li = local[it->second];
            addSum[li] += rec.sum;
            addCnt[li] += (int32_t)rec.count;
            st.seamContributions++;
          }
        }
        for (size_t i = 0; i < pts.size(); i++) {
          if (addCnt[i] == 0) continue;
          st.seamMatched++;
          vSum[pts[i]] += addSum[i];
          vCnt[pts[i]] += addCnt[i];
        }
      }

      void run(const std::string &caseDir, int procID)
      {
        load(caseDir, procID);

        // --- output vertices: all original points first
        out.vertices.resize(nP);
        for (int64_t i = 0; i < nP; i++)
          out.vertices[i] = umesh::vec3f((float)P[3*i], (float)P[3*i+1], (float)P[3*i+2]);
        out.values.assign(nP, 0.f);
        vSum.assign(nP, 0.0);
        vCnt.assign(nP, 0);
        faceCentreIdx.assign(nF, -1);

        const size_t expect = size_t(nC) * 11 / 10;
        out.hexes.reserve(expect);

        if (opt.wallMaskLayers > 0) computeWallMask();
        for (int64_t c = 0; c < nC; c++) convertCell(c);

        if (!opt.seamDir.empty()) applySeams(opt.seamDir);

        // --- cell -> vertex averaging (within this processor, plus the
        //     other processors' contributions if seams were applied)
        int64_t unused = 0;
        for (int64_t i = 0; i < nP; i++) {
          if (vCnt[i] > 0) out.values[i] = float(vSum[i] / vCnt[i]);
          else unused++;
        }
        for (auto &pf : pendingFaceCentres) {
          const int32_t f = pf.second;
          double s = 0.0;
          const int n = fOfs[f+1] - fOfs[f];
          for (int k = fOfs[f]; k < fOfs[f+1]; k++) s += out.values[fVtx[k]];
          out.values[pf.first] = float(s / n);
        }
        // wall mask: overwrite values of all vertices of masked cells
        for (size_t i = 0; i < vMask.size() && i < out.values.size(); i++)
          if (vMask[i]) { out.values[i] = opt.wallMaskValue; st.maskedVertices++; }
        if (unused)
          fprintf(stderr, "#f2u: processor%i: %lli unused points\n", procID, (long long)unused);

        st.outVertices = (int64_t)out.vertices.size();

        // --- bounds
        for (const auto &v : out.vertices) {
          const double p[3] = { v.x, v.y, v.z };
          out.bounds.extend(p);
        }

        // --- boundary patch bounds (aircraft framing)
        for (const auto &pa : patches) {
          if (pa.type == "processor" || pa.type == "processorCyclic") continue;
          out.patchTypes[pa.name] = pa.type;
          if (pa.nFaces <= 0) continue;
          auto &entry = out.patchBounds[pa.name];
          entry.first += pa.nFaces;
          for (int64_t f = pa.startFace; f < pa.startFace + pa.nFaces; f++)
            for (int k = fOfs[f]; k < fOfs[f+1]; k++)
              entry.second.extend(&P[3*size_t(fVtx[k])]);
        }

        out.vertices.shrink_to_fit();
        out.values.shrink_to_fit();
        out.tets.shrink_to_fit();
        out.pyrs.shrink_to_fit();
        out.wedges.shrink_to_fit();
        out.hexes.shrink_to_fit();
      }
    };

  } // anonymous

  void writeSeamFile(const std::string &caseDir, int procID,
                     const ConvertOptions &opt, ProcStats &stats)
  {
    ProcMesh tmp;
    Converter conv(opt, tmp);
    conv.load(caseDir, procID);
    conv.accumulateOnly();
    conv.writeSeam(opt.seamDir);
    stats = tmp.stats;
  }

  void convertProcessor(const std::string &caseDir, int procID,
                        const ConvertOptions &opt, ProcMesh &out)
  {
    out = ProcMesh();
    Converter conv(opt, out);
    conv.run(caseDir, procID);
  }

} // ::f2u
