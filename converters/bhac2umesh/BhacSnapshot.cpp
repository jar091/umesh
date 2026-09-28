// Reader for BHAC native snapshots; see BhacSnapshot.h.
// .par parsing, tail parsing and the forest walk are adapted from the
// SpaceConverter BHAC reader (space-converter/src/bhac/bhac_extract_iolib.cpp).

#include "BhacSnapshot.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace bhac {

  std::string lower(std::string s)
  {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return s;
  }

  // Fortran real literal ("0.5d0", "1.0D-3") -> double
  static double fortranDouble(std::string s)
  {
    for (char &c : s)
      if (c == 'd' || c == 'D') c = 'e';
    return std::stod(s);
  }

  // key = value pairs of a Fortran namelist file (keys lower case, quotes and
  // '!' comments removed). Values spread over several lines (e.g. typeB) are not
  // needed and only their first line is kept.
  static std::map<std::string, std::string> readNamelist(const std::string &path)
  {
    std::ifstream f(path);
    if (!f)
      throw std::runtime_error("BHAC: cannot open the .par file " + path);
    std::map<std::string, std::string> kv;
    const std::regex pairRe("([A-Za-z_][A-Za-z0-9_]*(?:\\([0-9, ]+\\))?)\\s*=\\s*('[^']*'|\"[^\"]*\"|[^,\\s]+)");
    std::string line;
    while (std::getline(f, line)) {
      bool quoted = false;
      char q = 0;
      for (size_t i = 0; i < line.size(); i++) {
        char c = line[i];
        if (quoted) {
          if (c == q) quoted = false;
        } else if (c == '\'' || c == '"') {
          quoted = true;
          q = c;
        } else if (c == '!') {
          line.erase(i);
          break;
        }
      }
      for (std::sregex_iterator it(line.begin(), line.end(), pairRe), end; it != end; ++it) {
        std::string key = lower((*it)[1].str());
        key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
        std::string val = (*it)[2].str();
        if (val.size() >= 2 && (val[0] == '\'' || val[0] == '"'))
          val = val.substr(1, val.size() - 2);
        kv[key] = val;
      }
    }
    return kv;
  }

  // Depth-first walk of the forest (forest.t read_node): leaves in file order
  static void walkNode(const std::vector<int32_t> &nodes, size_t &pos, int ndim, int level,
                       const int ig[3], std::vector<Leaf> &leaves)
  {
    if (pos >= nodes.size())
      throw std::runtime_error("BHAC: the forest ends early");
    const bool leaf = nodes[pos++] != 0;
    if (leaf) {
      leaves.push_back({ level, { ig[0], ig[1], ig[2] } });
      return;
    }
    const int n3 = ndim > 2 ? 2 : 1, n2 = ndim > 1 ? 2 : 1;
    for (int ic3 = 1; ic3 <= n3; ic3++)
      for (int ic2 = 1; ic2 <= n2; ic2++)
        for (int ic1 = 1; ic1 <= 2; ic1++) {
          const int ic[3] = { ic1, ic2, ic3 };
          int child[3];
          for (int d = 0; d < 3; d++)
            child[d] = (d < ndim) ? 2 * (ig[d] - 1) + ic[d] : 1;
          walkNode(nodes, pos, ndim, level + 1, child, leaves);
        }
  }

  void Snapshot::open(const std::string &_datFile, const std::string &parFile)
  {
    datFile = _datFile;
    std::ifstream f(datFile, std::ios::binary);
    if (!f)
      throw std::runtime_error("BHAC: cannot open " + datFile);
    f.seekg(0, std::ios::end);
    const size_t fileSize = (size_t)f.tellg();

    // ---- tail ----
    const size_t tailFixed = 8 * sizeof(int32_t) + sizeof(double);
    if (fileSize < tailFixed)
      throw std::runtime_error("BHAC: " + datFile + " is too short for a snapshot");
    int32_t v[8];
    f.seekg((std::streamoff)(fileSize - tailFixed));
    f.read(reinterpret_cast<char *>(v), sizeof(v));
    f.read(reinterpret_cast<char *>(&time), sizeof(double));
    nleafs = v[0]; levmax = v[1]; ndim = v[2]; ndir = v[3];
    nw = v[4]; nws = v[5]; neqpar = v[6]; iteration = v[7];
    if (!f || ndim < 1 || ndim > 3 || ndir < ndim || ndir > 3 || nw < 1 || nws < 0 ||
        neqpar < 0 || nleafs < 1 || levmax < 1 || levmax > 30)
      throw std::runtime_error("BHAC: " + datFile + " has no valid snapshot tail (only native-endian, "
                               "staggered-field .dat files are supported)");
    const size_t tailSize = tailFixed + (size_t)ndim * sizeof(int32_t) + (size_t)neqpar * sizeof(double);
    if (fileSize < tailSize)
      throw std::runtime_error("BHAC: " + datFile + " is too short for its header");
    f.seekg((std::streamoff)(fileSize - tailSize));
    int32_t n[3] = { 1, 1, 1 };
    f.read(reinterpret_cast<char *>(n), (std::streamsize)(ndim * sizeof(int32_t)));
    eqpar.assign((size_t)neqpar, 0.0);
    if (neqpar > 0)
      f.read(reinterpret_cast<char *>(eqpar.data()), (std::streamsize)(neqpar * sizeof(double)));
    if (!f)
      throw std::runtime_error("BHAC: cannot read the header of " + datFile);
    for (int d = 0; d < 3; d++)
      nx[d] = (d < ndim) ? n[d] : 1;

    // ---- .par ----
    std::map<std::string, std::string> kv = readNamelist(parFile);
    for (int d = 0; d < ndim; d++) {
      const std::string s = std::to_string(d + 1);
      if (!kv.count("nxlone" + s) || !kv.count("xprobmin" + s) || !kv.count("xprobmax" + s))
        throw std::runtime_error("BHAC: " + parFile + " lacks nxlone" + s + ", xprobmin" + s + " or xprobmax" + s);
      nxlone[d] = std::stoi(kv["nxlone" + s]);
      xprobmin[d] = fortranDouble(kv["xprobmin" + s]);
      xprobmax[d] = fortranDouble(kv["xprobmax" + s]);
    }
    typeaxial = kv.count("typeaxial") ? lower(kv["typeaxial"]) : std::string("slab");
    if (typeaxial == "spherical") {
      // amrio.t readparameters: xprob^LIM^DE = xprob^LIM^DE * two * dpi
      for (int d = 1; d < ndim; d++) {
        xprobmin[d] *= 2.0 * M_PI;
        xprobmax[d] *= 2.0 * M_PI;
      }
    }
    wnames.clear();
    if (kv.count("wnames")) {
      std::stringstream ss(kv["wnames"]);
      std::string w;
      while (ss >> w) wnames.push_back(w);
    }
    if ((int)wnames.size() != nw) {
      wnames.clear();
      for (int i = 0; i < nw; i++) wnames.push_back("w" + std::to_string(i + 1));
    }

    // ---- block layout ----
    size_t nstg = 1;
    ncellBlock = 1;
    for (int d = 0; d < ndim; d++) {
      ncellBlock *= (size_t)nx[d];
      nstg *= (size_t)nx[d] + 1;
    }
    blockBytes = ncellBlock * (size_t)nw * sizeof(double) + (nws > 0 ? nstg * (size_t)nws * sizeof(double) : 0);
    const size_t forestStart = blockBytes * (size_t)nleafs;
    if (forestStart + (size_t)nleafs * sizeof(int32_t) + tailSize > fileSize)
      throw std::runtime_error("BHAC: " + datFile + " is shorter than its " + std::to_string(nleafs) + " blocks");
    const size_t nnodes = (fileSize - tailSize - forestStart) / sizeof(int32_t);

    for (int d = 0; d < ndim; d++) {
      if (nxlone[d] % nx[d] != 0)
        throw std::runtime_error("BHAC: nxlone" + std::to_string(d + 1) + " is not a multiple of the block size");
      ng[d] = nxlone[d] / nx[d];
    }

    // ---- forest ----
    std::vector<int32_t> nodes(nnodes);
    f.seekg((std::streamoff)forestStart);
    f.read(reinterpret_cast<char *>(nodes.data()), (std::streamsize)(nnodes * sizeof(int32_t)));
    if (!f)
      throw std::runtime_error("BHAC: cannot read the forest of " + datFile);
    leaves.clear();
    size_t pos = 0;
    for (int ig3 = 1; ig3 <= ng[2]; ig3++)
      for (int ig2 = 1; ig2 <= ng[1]; ig2++)
        for (int ig1 = 1; ig1 <= ng[0]; ig1++) {
          const int ig[3] = { ig1, ig2, ig3 };
          walkNode(nodes, pos, ndim, 1, ig, leaves);
        }
    if (pos != nnodes || (int)leaves.size() != nleafs)
      throw std::runtime_error("BHAC: the forest of " + datFile + " does not match the base grid of the .par file (" +
                               std::to_string(leaves.size()) + " leaves of " + std::to_string(nleafs) + ")");
  }

  void Snapshot::readBlock(std::ifstream &f, size_t b, std::vector<double> &w) const
  {
    w.resize(ncellBlock * (size_t)nw);
    f.seekg((std::streamoff)(blockBytes * b));
    f.read(reinterpret_cast<char *>(w.data()), (std::streamsize)(w.size() * sizeof(double)));
    if (!f)
      throw std::runtime_error("BHAC: cannot read block " + std::to_string(b) + " of " + datFile);
  }

  int Snapshot::varIndex(const std::string &name) const
  {
    const std::string n = lower(name);
    for (int i = 0; i < (int)wnames.size(); i++)
      if (lower(wnames[i]) == n) return i;
    return -1;
  }

} // namespace bhac
