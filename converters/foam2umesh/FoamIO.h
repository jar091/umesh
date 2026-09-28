// foam2umesh - minimal reader for (decomposed) OpenFOAM polyMesh and
// volume fields written in *binary* format (label=32, scalar=64 or 32).
//
// Only what is needed for the conversion is implemented:
//   - FoamFile header (format, arch, class, note)
//   - binary List<T> / labelList / vectorField / faceCompactList
//   - polyBoundaryMesh ("boundary" file, ASCII dictionary)
//   - volScalarField / volVectorField "internalField nonuniform List<...>"
#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <map>

namespace foam {

  struct Header {
    std::string format;   // "binary" / "ascii"
    std::string arch;     // e.g. "LSB;label=32;scalar=64"
    std::string cls;      // e.g. "faceCompactList"
    std::string object;
    std::string note;
    int labelBytes  = 4;
    int scalarBytes = 8;
    /*! byte offset right behind the closing '}' of the FoamFile dict */
    size_t dataBegin = 0;
  };

  /*! whole file in memory */
  struct File {
    std::string       path;
    std::vector<char> data;
    Header            header;
  };

  /*! read a file into memory and parse its FoamFile header */
  void readFile(const std::string &path, File &file);

  /*! parse integers of the form "nPoints:123 nCells:456 ..." from the note */
  std::map<std::string,int64_t> parseNote(const std::string &note);

  /*! skip white space and C/C++ comments starting at 'pos' */
  size_t skipWs(const File &f, size_t pos);

  /*! read a binary list of 'elemBytes'-sized elements at 'pos'
      (form "<N>(<bytes>)" or "<N>{<one element>}"); returns the
      number of elements, fills 'out' with raw bytes, updates pos
      to behind the closing bracket. */
  size_t readBinaryListRaw(const File &f, size_t &pos, size_t elemBytes,
                           std::vector<char> &out);

  /*! labelList (32 bit labels only) */
  void readLabelList(const std::string &path, std::vector<int32_t> &out,
                     Header *hdr = nullptr);

  /*! vectorField (points) as flat xyz doubles */
  void readPoints(const std::string &path, std::vector<double> &xyz);

  /*! faceCompactList: offsets (nFaces+1) and flat vertex list */
  void readFaces(const std::string &path,
                 std::vector<int32_t> &offsets,
                 std::vector<int32_t> &verts);

  struct Patch {
    std::string name;
    std::string type;
    int64_t     nFaces    = 0;
    int64_t     startFace = 0;
    int         neighbProcNo = -1;  // processor patches only
    std::vector<std::string> inGroups;
  };

  /*! polyBoundaryMesh */
  void readBoundary(const std::string &path, std::vector<Patch> &patches);

  /*! internalField of a volScalarField (numComponents==1) or
      volVectorField (numComponents==3); values returned as doubles,
      numComponents per cell */
  void readInternalField(const std::string &path,
                         std::vector<double> &values,
                         int &numComponents);

} // ::foam
