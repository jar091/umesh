// foam2umesh - conversion of one OpenFOAM processor mesh into
// umesh-style elements (tets, pyramids, wedges, hexes) with
// per-vertex scalars.
//
// Vertex-order convention (umesh == VTK == what barney evaluates):
//
//   Tet  (v0,v1,v2,v3)        : dot(v3-v0, cross(v1-v0, v2-v0)) > 0,
//                               i.e. the right-hand normal of (v0,v1,v2)
//                               points towards v3. barney's tetScalar()
//                               rejects every point for tets with the
//                               opposite orientation (they would be
//                               invisible), so this one matters most.
//   Pyr  (v0..v3 | v4)        : right-hand normal of the base quad
//                               (v0,v1,v2,v3) points towards the apex v4.
//   Wedge(v0,v1,v2 | v3,v4,v5): right-hand normal of the triangle
//                               (v0,v1,v2) points AWAY from (v3,v4,v5)
//                               (VTK_WEDGE; see umesh FaceConn.cpp,
//                               fixNegativeVolumeElements.cpp,
//                               tetrahedralize.cpp); lateral edges are
//                               v0-v3, v1-v4, v2-v5.
//   Hex  (v0..v3 | v4..v7)    : right-hand normal of the base quad
//                               (v0,v1,v2,v3) points towards the top
//                               quad; vertical edges v0-v4 ... v3-v7.
//
// OpenFOAM faces: the right-hand normal of face f points out of
// owner[f] and into neighbour[f]; for the neighbour cell the face is
// therefore used in reversed order.
#pragma once

#include "umesh/UMesh.h"
#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <array>

namespace f2u {

  struct Box {
    double lo[3] = { +1e300, +1e300, +1e300 };
    double hi[3] = { -1e300, -1e300, -1e300 };
    bool empty() const { return lo[0] > hi[0]; }
    void extend(const double *p) {
      for (int k = 0; k < 3; k++) {
        lo[k] = std::min(lo[k], p[k]);
        hi[k] = std::max(hi[k], p[k]);
      }
    }
    void extend(const Box &b) {
      if (b.empty()) return;
      extend(b.lo); extend(b.hi);
    }
  };

  /*! statistics of one converted processor (all counters additive) */
  struct ProcStats {
    int64_t cells = 0;            // original OpenFOAM cells
    int64_t points = 0;           // original OpenFOAM points
    int64_t faces = 0;
    int64_t internalFaces = 0;
    int64_t facesBySize[3] = {0,0,0};  // 3, 4, >4 vertices

    // how original cells were emitted
    int64_t nativeHex = 0, nativeWedge = 0, nativePyr = 0, nativeTet = 0;
    int64_t decomposedPoly = 0;     // general polyhedra (split hexes, ...)
    int64_t decomposedQuality = 0;  // hex/wedge/pyr/tet shapes that were
                                    // decomposed because of a non-positive
                                    // corner Jacobian
    int64_t decomposedTopology = 0; // hex/wedge/pyr-like face sets whose
                                    // topology check failed

    // output
    int64_t outVertices = 0;
    int64_t addedCellCentres = 0;
    int64_t addedFaceCentres = 0;
    int64_t outTets = 0, outPyrs = 0, outWedges = 0, outHexes = 0;

    // checks of every emitted element (double precision, umesh
    // convention, see checkElement()):
    //  inverted*: signed volume <= 0, i.e. wrong orientation (must be 0)
    //  notStar* : at least one face-fan sub-tet towards the element
    //             centre is non-positive (folded/concave element)
    int64_t invertedTets = 0, invertedPyrs = 0, invertedWedges = 0, invertedHexes = 0;
    int64_t notStarTets = 0, notStarPyrs = 0, notStarWedges = 0, notStarHexes = 0;
    // native hexes/wedges with a non-positive corner Jacobian (kept or
    // decomposed depending on the quality mode)
    int64_t badCornerHexes = 0, badCornerWedges = 0;
    // quad faces of decomposed cells whose pyramid was not star-shaped
    // and that were split into two tets instead
    int64_t quadsSplitToTets = 0;

    // cross-processor vertex averaging (--seam-dir)
    int64_t seamPoints = 0;        // points on processor patches
    int64_t seamMatched = 0;       // ... that found >= 1 copy on another processor
    int64_t seamContributions = 0; // (point, other processor) pairs added

    // wall mask (--wall-mask)
    int64_t maskedCells = 0, maskedVertices = 0;

    int64_t nanCells = 0;
    double  cellMin = +1e300, cellMax = -1e300;

    void add(const ProcStats &o);
  };

  /*! the converted mesh of one OpenFOAM processor; vertex indices are
      local to this processor */
  struct ProcMesh {
    int procID = -1;
    std::vector<umesh::vec3f> vertices;
    std::vector<float>        values;
    std::vector<umesh::Tet>   tets;
    std::vector<umesh::Pyr>   pyrs;
    std::vector<umesh::Wedge> wedges;
    std::vector<umesh::Hex>   hexes;
    Box                       bounds;
    ProcStats                 stats;
    /*! bounds of the boundary faces of every non-processor patch
        with faces on this processor: name -> (nFaces, bounds) */
    std::map<std::string, std::pair<int64_t,Box>> patchBounds;
    std::map<std::string, std::string>            patchTypes;
    size_t memoryBytes() const;
  };

  struct ConvertOptions {
    std::string field = "p";
    std::string time  = "1";
    /*! when to decompose a native hex/wedge/pyr cell (towards its
        cell centre) instead of emitting it as such:
          0 = never,
          1 = if it is not star-shaped (see checkElement()),
          2 = if it is not star-shaped or has a non-positive corner
              Jacobian (hex/wedge) */
    int qualityMode = 1;
    /*! in decomposed cells, replace a non-star-shaped pyramid by two
        tets if that makes both tets positive */
    bool splitBadPyramids = true;
    /*! if not empty: directory with per-processor seam files written by
        writeSeamFile(); vertex values on processor boundaries then
        average the incident cells of ALL processors sharing the point
        (face neighbours and their face neighbours), so the duplicated
        boundary vertices of neighbouring processors/parts carry
        identical values (no seams). */
    std::string seamDir;
    /*! if > 0: all vertices of cells within this many cell layers of a
        wall patch get 'wallMaskValue' (layer 1 = cells owning a wall
        face, layer n+1 = face neighbours of layer n; per processor) */
    int   wallMaskLayers = 0;
    float wallMaskValue  = 0.f;
  };

  /*! seam pre-pass: for all points on processor patches of this
      processor, write (x,y,z, sum of incident cell values, number of
      incident cells) plus the list of neighbouring processors to
      <seamDir>/proc<procID>.seam */
  void writeSeamFile(const std::string &caseDir, int procID,
                     const ConvertOptions &opt, ProcStats &stats);

  /*! convert processor directory <caseDir>/processor<procID> */
  void convertProcessor(const std::string &caseDir, int procID,
                        const ConvertOptions &opt, ProcMesh &out);

  // ------------------------------------------------------------------
  // signed-volume checks (umesh/barney convention), shared by the
  // converter and the --check mode. Return true if all tested
  // sub-volumes are strictly positive.
  // ------------------------------------------------------------------
  template<typename A, typename B, typename C, typename D>
  inline double tetVol(const A &a, const B &b, const C &c, const D &d)
  {
    const double ux = double(b[0])-a[0], uy = double(b[1])-a[1], uz = double(b[2])-a[2];
    const double vx = double(c[0])-a[0], vy = double(c[1])-a[1], vz = double(c[2])-a[2];
    const double wx = double(d[0])-a[0], wy = double(d[1])-a[1], wz = double(d[2])-a[2];
    // dot(d-a, cross(b-a, c-a))
    return wx*(uy*vz-uz*vy) + wy*(uz*vx-ux*vz) + wz*(ux*vy-uy*vx);
  }

  template<typename V> bool tetOK(const V *v) { return tetVol(v[0],v[1],v[2],v[3]) > 0.0; }

  /*! corner-Jacobian test for hexes: base (v0..v3) normal points
      towards the top (v4..v7). Strict quality criterion. */
  template<typename V> bool hexCornersOK(const V *v)
  {
    static const int c[8][4] = {
      {0,1,3,4},{1,2,0,5},{2,3,1,6},{3,0,2,7},
      {4,7,5,0},{5,4,6,1},{6,5,7,2},{7,6,4,3} };
    for (int i = 0; i < 8; i++)
      if (!(tetVol(v[c[i][0]],v[c[i][1]],v[c[i][2]],v[c[i][3]]) > 0.0)) return false;
    return true;
  }

  /*! corner-Jacobian test for wedges ((v0,v1,v2) normal points away
      from (v3,v4,v5)) */
  template<typename V> bool wedgeCornersOK(const V *v)
  {
    static const int c[6][4] = {
      {0,2,1,3},{1,0,2,4},{2,1,0,5},{3,4,5,0},{4,5,3,1},{5,3,4,2} };
    for (int i = 0; i < 6; i++)
      if (!(tetVol(v[c[i][0]],v[c[i][1]],v[c[i][2]],v[c[i][3]]) > 0.0)) return false;
    return true;
  }

  /*! INWARD facing faces of each element type in the umesh convention
      (identical to umesh/FaceConn.cpp); -1 terminates triangles */
  struct ElementFaces { int numVerts, numFaces; int f[6][4]; };
  static const ElementFaces tetFaces   = { 4, 4, {{1,3,2,-1},{0,2,3,-1},{0,3,1,-1},{0,1,2,-1}} };
  static const ElementFaces pyrFaces   = { 5, 5, {{4,1,0,-1},{4,2,1,-1},{4,3,2,-1},{4,0,3,-1},{0,1,2,3}} };
  static const ElementFaces wedgeFaces = { 6, 5, {{0,2,1,-1},{3,4,5,-1},{0,3,5,2},{1,2,5,4},{0,1,4,3}} };
  static const ElementFaces hexFaces   = { 8, 6, {{0,1,2,3},{4,7,6,5},{0,4,5,1},{2,6,7,3},{1,5,6,2},{0,3,7,4}} };

  /*! orientation / star-shape test of an element in the umesh
      convention: every face (inward orientation) is fanned from its
      vertex average fc and each triangle (fc,a,b) forms a sub-tet with
      the element's vertex average c. 'total' is the signed volume of
      the element: > 0 iff the vertex order has the correct
      orientation. 'allPositive' additionally requires every sub-tet to
      be positive (element star-shaped w.r.t. its centre = no folded
      faces). */
  struct VolCheck { double total = 0.0; bool allPositive = true; };

  template<typename V>
  inline VolCheck checkElement(const V *v, const ElementFaces &ef)
  {
    VolCheck r;
    double c[3] = {0,0,0};
    for (int i = 0; i < ef.numVerts; i++) for (int k = 0; k < 3; k++) c[k] += v[i][k];
    for (int k = 0; k < 3; k++) c[k] /= ef.numVerts;
    for (int fi = 0; fi < ef.numFaces; fi++) {
      const int *f = ef.f[fi];
      const int n = (f[3] < 0) ? 3 : 4;
      if (n == 3) {
        const double vol = tetVol(v[f[0]], v[f[1]], v[f[2]], c);
        r.total += vol;
        if (!(vol > 0.0)) r.allPositive = false;
        continue;
      }
      double fc[3] = {0,0,0};
      for (int i = 0; i < n; i++) for (int k = 0; k < 3; k++) fc[k] += v[f[i]][k];
      for (int k = 0; k < 3; k++) fc[k] /= n;
      for (int i = 0; i < n; i++) {
        const double vol = tetVol(fc, v[f[i]], v[f[(i+1)%n]], c);
        r.total += vol;
        if (!(vol > 0.0)) r.allPositive = false;
      }
    }
    return r;
  }

} // ::f2u
