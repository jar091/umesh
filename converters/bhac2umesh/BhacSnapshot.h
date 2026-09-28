// Reader for BHAC (MPI-AMRVAC family) native snapshots (dataNNNN.dat).
//
// The file layout, the .par parsing and the forest walk are adapted from the
// SpaceConverter BHAC reader (space-converter/src/bhac/bhac_extract_iolib.cpp,
// IT4Innovations, GPL-3.0-or-later).
//
// File format (BHAC amrio.t write_snapshot, native endianness, no record markers):
//   nleafs x { w(nx1,nx2,nx3,1:nw) [+ ws(0:nx1,0:nx2,0:nx3,1:nws) if staggered] }
//       blocks in Morton order = the order of the leaves in the forest below,
//       each variable a Fortran-ordered (x1 fastest) array of doubles without
//       ghost cells (MPI subarray type_block_io = ixM range, all nw variables)
//   forest: one 4-byte logical per tree node, depth first; roots in ig1-fastest
//       order over ng1 x ng2 x ng3 level-1 blocks, children in ic1-fastest order
//   nx(1:ndim) int32, eqpar(1:neqpar+nspecialpar) double,
//   nleafs, levmax, ndim, ndir, nw, nws, neqpar+nspecialpar, it  int32, t double
// The base grid and the domain are not in the file: they come from the .par file
// (nxlone^D, xprobmin^D, xprobmax^D, typeaxial; spherical: x2, x3 in units of 2 pi).
#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace bhac {

  struct Leaf {
    int level;
    int ig[3];   // 1-based block index on its level
  };

  struct Snapshot {
    // .par settings
    int nxlone[3] = { 1, 1, 1 };
    double xprobmin[3] = { 0.0, 0.0, 0.0 };
    double xprobmax[3] = { 1.0, 1.0, 1.0 };   // spherical: x2, x3 already multiplied by 2 pi
    std::string typeaxial = "slab";
    std::vector<std::string> wnames;

    // .dat tail
    int ndim = 3, ndir = 3, nw = 0, nws = 0, neqpar = 0, nleafs = 0, levmax = 1, iteration = 0;
    int nx[3] = { 1, 1, 1 };
    double time = 0.0;
    std::vector<double> eqpar;

    // derived layout
    int ng[3] = { 1, 1, 1 };            // base grid in blocks
    size_t ncellBlock = 0;              // cells per block
    size_t blockBytes = 0;              // bytes per block incl. staggered field
    std::vector<Leaf> leaves;           // file order
    std::string datFile;

    /*! reads the tail, the forest and the .par file; throws on errors */
    void open(const std::string &datFile, const std::string &parFile);

    /*! reads the cell-centred variables of leaf block b: w[v*ncellBlock + c],
        c = i + nx1*(j + nx2*k) */
    void readBlock(std::ifstream &f, size_t b, std::vector<double> &w) const;

    /*! index of a variable by (case insensitive) name, -1 if not present */
    int varIndex(const std::string &name) const;
  };

  std::string lower(std::string s);

} // namespace bhac
