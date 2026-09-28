// foam2umesh - minimal binary OpenFOAM reader (see FoamIO.h)

#include "FoamIO.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <fstream>

namespace foam {

  static std::runtime_error error(const std::string &path, const std::string &msg)
  {
    return std::runtime_error("foam: " + path + ": " + msg);
  }

  // ------------------------------------------------------------------
  // generic helpers
  // ------------------------------------------------------------------

  size_t skipWs(const File &f, size_t pos)
  {
    const std::vector<char> &d = f.data;
    const size_t n = d.size();
    while (pos < n) {
      const char c = d[pos];
      if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v') {
        ++pos;
      } else if (c == '/' && pos + 1 < n && d[pos+1] == '/') {
        while (pos < n && d[pos] != '\n') ++pos;
      } else if (c == '/' && pos + 1 < n && d[pos+1] == '*') {
        pos += 2;
        while (pos + 1 < n && !(d[pos] == '*' && d[pos+1] == '/')) ++pos;
        pos += 2;
      } else {
        break;
      }
    }
    return pos;
  }

  /*! read a word token (anything up to white space or one of ";{}()[]") */
  static std::string readWord(const File &f, size_t &pos)
  {
    pos = skipWs(f, pos);
    const std::vector<char> &d = f.data;
    std::string w;
    if (pos < d.size() && d[pos] == '"') {
      ++pos;
      while (pos < d.size() && d[pos] != '"') w.push_back(d[pos++]);
      ++pos;
      return w;
    }
    while (pos < d.size()) {
      const char c = d[pos];
      if (isspace((unsigned char)c) || strchr(";{}()[]", c)) break;
      w.push_back(c);
      ++pos;
    }
    return w;
  }

  static int64_t readInt(const File &f, size_t &pos)
  {
    const std::string w = readWord(f, pos);
    if (w.empty())
      throw error(f.path, "expected integer at byte " + std::to_string(pos));
    char *end = nullptr;
    const long long v = strtoll(w.c_str(), &end, 10);
    if (*end != 0)
      throw error(f.path, "expected integer, got '" + w + "'");
    return v;
  }

  // ------------------------------------------------------------------
  // header
  // ------------------------------------------------------------------

  static void parseHeader(File &f)
  {
    const std::vector<char> &d = f.data;
    const char *key = "FoamFile";
    const size_t searchLen = std::min<size_t>(d.size(), 1 << 16);
    const std::string head(d.data(), searchLen);
    size_t p = head.find(key);
    if (p == std::string::npos)
      throw error(f.path, "no FoamFile header");
    p = head.find('{', p);
    const size_t e = head.find('}', p);
    if (p == std::string::npos || e == std::string::npos)
      throw error(f.path, "broken FoamFile header");
    // tokenise "key value;" entries
    std::string body = head.substr(p + 1, e - p - 1);
    std::istringstream in(body);
    std::string line;
    // entries are separated by ';'
    size_t s = 0;
    while (s < body.size()) {
      size_t semi = body.find(';', s);
      if (semi == std::string::npos) break;
      std::string entry = body.substr(s, semi - s);
      s = semi + 1;
      // trim
      size_t a = entry.find_first_not_of(" \t\r\n");
      if (a == std::string::npos) continue;
      entry = entry.substr(a);
      size_t b = entry.find_first_of(" \t\r\n");
      if (b == std::string::npos) continue;
      std::string k = entry.substr(0, b);
      std::string v = entry.substr(b);
      a = v.find_first_not_of(" \t\r\n");
      v = (a == std::string::npos) ? "" : v.substr(a);
      while (!v.empty() && isspace((unsigned char)v.back())) v.pop_back();
      if (v.size() >= 2 && v.front() == '"' && v.back() == '"')
        v = v.substr(1, v.size() - 2);
      if (k == "format") f.header.format = v;
      else if (k == "arch") f.header.arch = v;
      else if (k == "class") f.header.cls = v;
      else if (k == "object") f.header.object = v;
      else if (k == "note") f.header.note = v;
    }
    f.header.dataBegin = e + 1;

    // arch: "LSB;label=32;scalar=64"
    const std::string &arch = f.header.arch;
    if (!arch.empty()) {
      if (arch.find("MSB") != std::string::npos)
        throw error(f.path, "big-endian files not supported");
      size_t l = arch.find("label=");
      if (l != std::string::npos) f.header.labelBytes = atoi(arch.c_str() + l + 6) / 8;
      size_t sc = arch.find("scalar=");
      if (sc != std::string::npos) f.header.scalarBytes = atoi(arch.c_str() + sc + 7) / 8;
    }
    if (f.header.format != "binary")
      throw error(f.path, "only binary format is supported (got '" + f.header.format + "')");
    if (f.header.labelBytes != 4)
      throw error(f.path, "only label=32 is supported");
    if (f.header.scalarBytes != 8 && f.header.scalarBytes != 4)
      throw error(f.path, "unsupported scalar size");
  }

  void readFile(const std::string &path, File &file)
  {
    file.path = path;
    FILE *fp = fopen(path.c_str(), "rb");
    if (!fp) throw error(path, "cannot open file");
    fseek(fp, 0, SEEK_END);
    const long size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    file.data.resize(size);
    if (size > 0 && fread(file.data.data(), 1, size, fp) != (size_t)size) {
      fclose(fp);
      throw error(path, "short read");
    }
    fclose(fp);
    parseHeader(file);
  }

  std::map<std::string,int64_t> parseNote(const std::string &note)
  {
    std::map<std::string,int64_t> r;
    std::istringstream in(note);
    std::string tok;
    while (in >> tok) {
      size_t c = tok.find(':');
      if (c == std::string::npos) continue;
      r[tok.substr(0, c)] = atoll(tok.c_str() + c + 1);
    }
    return r;
  }

  // ------------------------------------------------------------------
  // binary lists
  // ------------------------------------------------------------------

  size_t readBinaryListRaw(const File &f, size_t &pos, size_t elemBytes,
                           std::vector<char> &out)
  {
    const int64_t n = readInt(f, pos);
    if (n < 0) throw error(f.path, "negative list size");
    pos = skipWs(f, pos);
    if (pos >= f.data.size()) throw error(f.path, "unexpected EOF");
    const char open = f.data[pos];
    if (open == '(') {
      ++pos;
      const size_t bytes = size_t(n) * elemBytes;
      if (pos + bytes > f.data.size()) throw error(f.path, "list exceeds file size");
      out.assign(f.data.begin() + pos, f.data.begin() + pos + bytes);
      pos += bytes;
      if (pos >= f.data.size() || f.data[pos] != ')')
        throw error(f.path, "expected ')' after binary list");
      ++pos;
    } else if (open == '{') {
      // uniform list: one element, repeated n times
      ++pos;
      if (pos + elemBytes > f.data.size()) throw error(f.path, "unexpected EOF");
      out.resize(size_t(n) * elemBytes);
      for (int64_t i = 0; i < n; i++)
        memcpy(out.data() + i * elemBytes, f.data.data() + pos, elemBytes);
      pos += elemBytes;
      pos = skipWs(f, pos);
      if (pos >= f.data.size() || f.data[pos] != '}')
        throw error(f.path, "expected '}' after uniform list");
      ++pos;
    } else {
      throw error(f.path, std::string("expected '(' or '{' after list size, got '") + open + "'");
    }
    return size_t(n);
  }

  template<typename T>
  static void rawTo(const std::vector<char> &raw, std::vector<T> &out)
  {
    out.resize(raw.size() / sizeof(T));
    if (!raw.empty()) memcpy(out.data(), raw.data(), out.size() * sizeof(T));
  }

  static void rawScalarsToDouble(const std::vector<char> &raw, int scalarBytes,
                                 std::vector<double> &out)
  {
    if (scalarBytes == 8) { rawTo(raw, out); return; }
    std::vector<float> tmp;
    rawTo(raw, tmp);
    out.assign(tmp.begin(), tmp.end());
  }

  void readLabelList(const std::string &path, std::vector<int32_t> &out, Header *hdr)
  {
    File f;
    readFile(path, f);
    size_t pos = f.header.dataBegin;
    std::vector<char> raw;
    readBinaryListRaw(f, pos, 4, raw);
    rawTo(raw, out);
    if (hdr) *hdr = f.header;
  }

  void readPoints(const std::string &path, std::vector<double> &xyz)
  {
    File f;
    readFile(path, f);
    size_t pos = f.header.dataBegin;
    std::vector<char> raw;
    readBinaryListRaw(f, pos, 3 * f.header.scalarBytes, raw);
    rawScalarsToDouble(raw, f.header.scalarBytes, xyz);
  }

  void readFaces(const std::string &path,
                 std::vector<int32_t> &offsets,
                 std::vector<int32_t> &verts)
  {
    File f;
    readFile(path, f);
    if (f.header.cls != "faceCompactList")
      throw error(path, "expected class faceCompactList, got '" + f.header.cls + "'");
    size_t pos = f.header.dataBegin;
    std::vector<char> raw;
    readBinaryListRaw(f, pos, 4, raw);
    rawTo(raw, offsets);
    readBinaryListRaw(f, pos, 4, raw);
    rawTo(raw, verts);
    if (offsets.empty() || (size_t)offsets.back() != verts.size())
      throw error(path, "inconsistent faceCompactList");
  }

  // ------------------------------------------------------------------
  // boundary
  // ------------------------------------------------------------------

  void readBoundary(const std::string &path, std::vector<Patch> &patches)
  {
    File f;
    readFile(path, f);
    size_t pos = f.header.dataBegin;
    const int64_t n = readInt(f, pos);
    pos = skipWs(f, pos);
    if (f.data[pos] != '(') throw error(path, "expected '('");
    ++pos;
    patches.clear();
    for (int64_t i = 0; i < n; i++) {
      Patch p;
      p.name = readWord(f, pos);
      pos = skipWs(f, pos);
      if (f.data[pos] != '{') throw error(path, "expected '{' after patch name");
      ++pos;
      while (true) {
        pos = skipWs(f, pos);
        if (f.data[pos] == '}') { ++pos; break; }
        std::string key = readWord(f, pos);
        // value: everything up to ';' (may contain "1(wall)")
        std::string val;
        while (pos < f.data.size() && f.data[pos] != ';') val.push_back(f.data[pos++]);
        ++pos;
        size_t a = val.find_first_not_of(" \t\r\n");
        val = (a == std::string::npos) ? "" : val.substr(a);
        while (!val.empty() && isspace((unsigned char)val.back())) val.pop_back();
        if (key == "type") p.type = val;
        else if (key == "nFaces") p.nFaces = atoll(val.c_str());
        else if (key == "startFace") p.startFace = atoll(val.c_str());
        else if (key == "neighbProcNo") p.neighbProcNo = atoi(val.c_str());
        else if (key == "inGroups") {
          size_t o = val.find('('), c = val.find(')');
          if (o != std::string::npos && c != std::string::npos) {
            std::istringstream gs(val.substr(o + 1, c - o - 1));
            std::string g;
            while (gs >> g) p.inGroups.push_back(g);
          }
        }
      }
      patches.push_back(p);
    }
  }

  // ------------------------------------------------------------------
  // fields
  // ------------------------------------------------------------------

  void readInternalField(const std::string &path,
                         std::vector<double> &values,
                         int &numComponents)
  {
    File f;
    readFile(path, f);
    // the text "internalField" appears before any binary payload, so a
    // plain search is safe
    const char *key = "internalField";
    const std::vector<char> &d = f.data;
    size_t pos = std::string::npos;
    for (size_t i = f.header.dataBegin; i + strlen(key) < d.size(); i++)
      if (memcmp(d.data() + i, key, strlen(key)) == 0) { pos = i + strlen(key); break; }
    if (pos == std::string::npos) throw error(path, "no internalField");

    const std::string kind = readWord(f, pos);
    if (kind != "nonuniform")
      throw error(path, "internalField is '" + kind +
                  "' - only 'nonuniform List<...>' is supported (use a written time step)");
    const std::string type = readWord(f, pos);
    if (type == "List<scalar>") numComponents = 1;
    else if (type == "List<vector>") numComponents = 3;
    else throw error(path, "unsupported internalField type '" + type + "'");

    std::vector<char> raw;
    readBinaryListRaw(f, pos, numComponents * f.header.scalarBytes, raw);
    rawScalarsToDouble(raw, f.header.scalarBytes, values);
  }

} // ::foam
