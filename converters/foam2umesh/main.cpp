// foam2umesh - convert a decomposed OpenFOAM case (binary polyMesh in
// processorN/constant/polyMesh, cell-centred volScalar/VectorFields)
// into N stand-alone .umesh part files for data-parallel rendering with
// haystack/barney.
//
// Modes:
//   convert   (default) : converts one or more "blocks" of processors and
//                         writes all requested part files that fall into
//                         those blocks, plus per-block statistics.
//   --summarize         : merges the per-block statistics into a JSON/txt
//                         summary (bounds, value range, percentiles, part
//                         sizes, element totals, GPU memory estimate).
//   --check <f.umesh>   : reloads a part with UMesh::loadFrom and validates it.

#include "FoamIO.h"
#include "Convert.h"

#include <omp.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <array>
#include <initializer_list>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace f2u;

namespace {

  double now()
  {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
  }

  void usage(const std::string &err = "")
  {
    if (!err.empty()) std::cerr << "error: " << err << "\n\n";
    std::cerr <<
      "usage:\n"
      "  foam2umesh --case <dir> -n <N>[,<N2>,...] -o <base> [options]\n"
      "      converts processor0..P-1 of a decomposed OpenFOAM case into\n"
      "      N parts; part i is written to <base>_<iiii>.umesh. <base> may\n"
      "      contain '{N}' (required when several N are given).\n"
      "    --field <name>      field in <time>/ (default p); vector fields -> magnitude\n"
      "    --time <dir>|latest time directory (default latest)\n"
      "    --nprocs <P>        number of processor dirs (default: auto)\n"
      "    --stats-dir <dir>   per-block statistics (default <dir of base>/stats)\n"
      "    --grouping contiguous|hierarchical:nx,ny,nz\n"
      "                        how processors are grouped into parts: contiguous\n"
      "                        processor-id ranges (default), or subtrees of an\n"
      "                        OpenFOAM hierarchical (nx ny nz) decomposition\n"
      "                        (compact, nearly disjoint part boxes)\n"
      "    --blocks-begin <b> --blocks-end <e>\n"
      "                        range of blocks (block = 1024/min(N) processors)\n"
      "                        handled by this srun step; tasks of the step\n"
      "                        (SLURM_PROCID/SLURM_NTASKS) share them round-robin\n"
      "                        (default: all blocks)\n"
      "    --threads <t>       OpenMP threads (default: OMP_NUM_THREADS / all)\n"
      "    --quality <0|1|2>   decompose native hex/wedge/pyr cells to their centre\n"
      "                        0: never, 1: if not star-shaped (default),\n"
      "                        2: also if a hex/wedge corner Jacobian is <= 0\n"
      "    --no-split-pyramids keep non-star-shaped pyramids of decomposed cells\n"
      "                        (default: split them into 2 tets when that helps)\n"
      "    --sample-stride <k> keep every k-th vertex value for percentiles (128)\n"
      "    --seam-dir <dir>    average vertex values across processor boundaries using\n"
      "                        the seam files in <dir> (written by --seam-prepass)\n"
      "    --seam-prepass      only write the seam files (run before the conversion)\n"
      "    --wall-mask <k>     set the values of all vertices of cells within k cell\n"
      "                        layers of a wall patch (per processor) to\n"
      "    --wall-mask-value <v>  (default 0)\n"
      "  foam2umesh --summarize --stats-dir <dir> --summary <prefix> [--region name:x0,x1,y0,y1,z0,z1 ...]\n"
      "  foam2umesh --check <file.umesh>\n";
    exit(err.empty() ? 0 : 1);
  }

  std::vector<int> parseIntList(const std::string &s)
  {
    std::vector<int> r;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, ',')) if (!t.empty()) r.push_back(std::stoi(t));
    return r;
  }

  bool isDir(const std::string &p)
  {
    struct stat st;
    return stat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
  }

  void mkdirs(const std::string &path)
  {
    std::string cur;
    std::stringstream ss(path);
    std::string part;
    if (!path.empty() && path[0] == '/') cur = "/";
    while (std::getline(ss, part, '/')) {
      if (part.empty()) continue;
      cur += part + "/";
      mkdir(cur.c_str(), 0775);
    }
  }

  std::string dirName(const std::string &p)
  {
    const size_t s = p.rfind('/');
    return s == std::string::npos ? "." : p.substr(0, s);
  }

  std::string replaceAll(std::string s, const std::string &from, const std::string &to)
  {
    size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
      s.replace(pos, from.size(), to);
      pos += to.size();
    }
    return s;
  }

  std::string partFileName(const std::string &base, int N, int i)
  {
    char suffix[32];
    snprintf(suffix, sizeof(suffix), "_%04d.umesh", i);
    return replaceAll(base, "{N}", std::to_string(N)) + suffix;
  }

  int countProcessors(const std::string &caseDir)
  {
    int n = 0;
    while (isDir(caseDir + "/processor" + std::to_string(n))) n++;
    return n;
  }

  std::string latestTime(const std::string &procDir)
  {
    DIR *d = opendir(procDir.c_str());
    if (!d) throw std::runtime_error("cannot list " + procDir);
    std::string best;
    double bestT = -1e300;
    while (dirent *e = readdir(d)) {
      const std::string name = e->d_name;
      char *end = nullptr;
      const double t = strtod(name.c_str(), &end);
      if (name.empty() || *end != 0) continue;
      if (!isDir(procDir + "/" + name)) continue;
      if (t > bestT) { bestT = t; best = name; }
    }
    closedir(d);
    if (best.empty()) throw std::runtime_error("no time directories in " + procDir);
    return best;
  }

  // ------------------------------------------------------------------
  // tiny JSON helpers
  // ------------------------------------------------------------------
  std::string jnum(double v)
  {
    if (!std::isfinite(v)) return "null";
    char buf[64];
    snprintf(buf, sizeof(buf), "%.9g", v);
    return buf;
  }
  std::string jstr(const std::string &s) { return "\"" + s + "\""; }
  std::string jbox(const Box &b)
  {
    return "{\"lower\": [" + jnum(b.lo[0]) + ", " + jnum(b.lo[1]) + ", " + jnum(b.lo[2]) +
      "], \"upper\": [" + jnum(b.hi[0]) + ", " + jnum(b.hi[1]) + ", " + jnum(b.hi[2]) + "]}";
  }

  // ------------------------------------------------------------------
  // per-part record (also stored in the block stats files)
  // ------------------------------------------------------------------
  struct PartRecord {
    int N = 0, index = 0;
    std::vector<int> procs;   // OpenFOAM processor ids in this part
    int64_t cells = 0, vertices = 0, tets = 0, pyrs = 0, wedges = 0, hexes = 0;
    int64_t indices = 0, fileBytes = 0;
    Box bounds;
    double vmin = +1e300, vmax = -1e300;
    std::string file;
  };

  void writeBox(std::ostream &o, const Box &b)
  {
    o << b.lo[0] << ' ' << b.lo[1] << ' ' << b.lo[2] << ' '
      << b.hi[0] << ' ' << b.hi[1] << ' ' << b.hi[2];
  }
  void readBox(std::istream &i, Box &b)
  {
    i >> b.lo[0] >> b.lo[1] >> b.lo[2] >> b.hi[0] >> b.hi[1] >> b.hi[2];
  }

  void writeStats(std::ostream &o, const ProcStats &s)
  {
    o << s.cells << ' ' << s.points << ' ' << s.faces << ' ' << s.internalFaces << ' '
      << s.facesBySize[0] << ' ' << s.facesBySize[1] << ' ' << s.facesBySize[2] << ' '
      << s.nativeHex << ' ' << s.nativeWedge << ' ' << s.nativePyr << ' ' << s.nativeTet << ' '
      << s.decomposedPoly << ' ' << s.decomposedQuality << ' ' << s.decomposedTopology << ' '
      << s.outVertices << ' ' << s.addedCellCentres << ' ' << s.addedFaceCentres << ' '
      << s.outTets << ' ' << s.outPyrs << ' ' << s.outWedges << ' ' << s.outHexes << ' '
      << s.invertedTets << ' ' << s.invertedPyrs << ' ' << s.invertedWedges << ' ' << s.invertedHexes << ' '
      << s.notStarTets << ' ' << s.notStarPyrs << ' ' << s.notStarWedges << ' ' << s.notStarHexes << ' '
      << s.badCornerHexes << ' ' << s.badCornerWedges << ' ' << s.quadsSplitToTets << ' '
      << s.seamPoints << ' ' << s.seamMatched << ' ' << s.seamContributions << ' '
      << s.maskedCells << ' ' << s.maskedVertices << ' '
      << s.nanCells << ' ' << s.cellMin << ' ' << s.cellMax;
  }
  void readStats(std::istream &i, ProcStats &s)
  {
    i >> s.cells >> s.points >> s.faces >> s.internalFaces
      >> s.facesBySize[0] >> s.facesBySize[1] >> s.facesBySize[2]
      >> s.nativeHex >> s.nativeWedge >> s.nativePyr >> s.nativeTet
      >> s.decomposedPoly >> s.decomposedQuality >> s.decomposedTopology
      >> s.outVertices >> s.addedCellCentres >> s.addedFaceCentres
      >> s.outTets >> s.outPyrs >> s.outWedges >> s.outHexes
      >> s.invertedTets >> s.invertedPyrs >> s.invertedWedges >> s.invertedHexes
      >> s.notStarTets >> s.notStarPyrs >> s.notStarWedges >> s.notStarHexes
      >> s.badCornerHexes >> s.badCornerWedges >> s.quadsSplitToTets
      >> s.seamPoints >> s.seamMatched >> s.seamContributions
      >> s.maskedCells >> s.maskedVertices
      >> s.nanCells >> s.cellMin >> s.cellMax;
  }

  int64_t fileSize(const std::string &p)
  {
    struct stat st;
    return stat(p.c_str(), &st) == 0 ? (int64_t)st.st_size : -1;
  }

  // ------------------------------------------------------------------
  // assemble one part from converted processors and write it
  // ------------------------------------------------------------------
  PartRecord writePart(const std::vector<const ProcMesh *> &procs,
                       const std::string &fileName, const std::string &attrName)
  {
    PartRecord rec;
    umesh::UMesh mesh;
    size_t nv = 0, nt = 0, np = 0, nw = 0, nh = 0;
    for (auto p : procs) {
      nv += p->vertices.size(); nt += p->tets.size(); np += p->pyrs.size();
      nw += p->wedges.size(); nh += p->hexes.size();
      rec.cells += p->stats.cells;
    }
    if (nv >= (size_t)std::numeric_limits<int32_t>::max())
      throw std::runtime_error("part too large for 32-bit vertex indices: " + fileName);
    mesh.vertices.reserve(nv);
    mesh.tets.reserve(nt); mesh.pyrs.reserve(np); mesh.wedges.reserve(nw); mesh.hexes.reserve(nh);
    auto attr = std::make_shared<umesh::Attribute>();
    attr->name = attrName;
    attr->values.reserve(nv);

    for (auto p : procs) {
      const int32_t o = (int32_t)mesh.vertices.size();
      mesh.vertices.insert(mesh.vertices.end(), p->vertices.begin(), p->vertices.end());
      attr->values.insert(attr->values.end(), p->values.begin(), p->values.end());
      for (auto e : p->tets)   { for (int k = 0; k < 4; k++) e[k] += o; mesh.tets.push_back(e); }
      for (auto e : p->pyrs)   { for (int k = 0; k < 5; k++) e[k] += o; mesh.pyrs.push_back(e); }
      for (auto e : p->wedges) { for (int k = 0; k < 6; k++) e[k] += o; mesh.wedges.push_back(e); }
      for (auto e : p->hexes)  { for (int k = 0; k < 8; k++) e[k] += o; mesh.hexes.push_back(e); }
      rec.bounds.extend(p->bounds);
    }
    // bounds and value range directly (all vertices are referenced by
    // elements, so this equals what UMesh::finalize() would compute)
    umesh::box3f bb;
    bb.lower = umesh::vec3f((float)rec.bounds.lo[0], (float)rec.bounds.lo[1], (float)rec.bounds.lo[2]);
    bb.upper = umesh::vec3f((float)rec.bounds.hi[0], (float)rec.bounds.hi[1], (float)rec.bounds.hi[2]);
    mesh.bounds = bb;
    umesh::range1f vr;
    for (float v : attr->values) { vr.lower = std::min(vr.lower, v); vr.upper = std::max(vr.upper, v); }
    attr->valueRange = vr;
    mesh.perVertex = attr;

    const std::string tmp = fileName + ".tmp";
    mesh.saveTo(tmp);
    if (rename(tmp.c_str(), fileName.c_str()) != 0)
      throw std::runtime_error("cannot rename " + tmp);

    rec.vertices = nv; rec.tets = nt; rec.pyrs = np; rec.wedges = nw; rec.hexes = nh;
    rec.indices = 4*nt + 5*np + 6*nw + 8*nh;
    rec.vmin = vr.lower; rec.vmax = vr.upper;
    rec.file = fileName;
    rec.fileBytes = fileSize(fileName);
    return rec;
  }

  // ==================================================================
  // convert mode
  // ==================================================================
  struct Args {
    std::string caseDir, base, statsDir, summary, check, seamDir;
    bool seamPrepass = false;
    std::vector<std::string> regions;
    int wallMask = 0;
    float wallMaskValue = 0.f;
    std::string grouping = "contiguous";
    std::string field = "p", time = "latest";
    std::vector<int> Ns;
    int nprocs = 0, threads = 0;
    int blocksBegin = -1, blocksEnd = -1;
    int sampleStride = 128;
    int qualityMode = 1;
    bool splitBadPyramids = true;
    bool summarize = false;
  };

  int runConvert(Args &a)
  {
    if (a.caseDir.empty()) usage("--case missing");
    if (a.Ns.empty()) usage("-n missing");
    if (a.base.empty()) usage("-o missing");
    if (a.Ns.size() > 1 && a.base.find("{N}") == std::string::npos)
      usage("-o must contain {N} when several part counts are requested");
    if (a.threads > 0) omp_set_num_threads(a.threads);

    const int P = a.nprocs > 0 ? a.nprocs : countProcessors(a.caseDir);
    if (P <= 0) throw std::runtime_error("no processor directories in " + a.caseDir);
    const int minN = *std::min_element(a.Ns.begin(), a.Ns.end());
    for (int N : a.Ns)
      if (P % N != 0) throw std::runtime_error("number of processors " + std::to_string(P)
                                               + " not divisible by N=" + std::to_string(N));
    // processor order; parts are contiguous ranges of this order
    std::vector<int> order(P);
    for (int i = 0; i < P; i++) order[i] = i;
    if (a.grouping.rfind("hierarchical:", 0) == 0) {
      // OpenFOAM 'hierarchical' decomposition with n = (nx ny nz): the
      // processor id is xi + nx*(yi + ny*zi); the domain is cut into nx
      // x-slabs, each slab into ny y-columns, each column into nz
      // z-boxes. Consecutive z-boxes of one column (and neighbouring
      // columns of one slab) form disjoint boxes, so parts are taken as
      // contiguous ranges of the order (xi, yi, zi) with zi fastest.
      std::vector<int> n = parseIntList(a.grouping.substr(13));
      if (n.size() != 3 || n[0]*n[1]*n[2] != P)
        throw std::runtime_error("--grouping hierarchical:nx,ny,nz must multiply to " + std::to_string(P));
      int k = 0;
      for (int xi = 0; xi < n[0]; xi++)
        for (int yi = 0; yi < n[1]; yi++)
          for (int zi = 0; zi < n[2]; zi++)
            order[k++] = xi + n[0]*(yi + n[1]*zi);
    } else if (a.grouping != "contiguous") {
      throw std::runtime_error("unknown --grouping '" + a.grouping + "'");
    }
    const int B = P / minN;           // processors per block
    const int numBlocks = minN;
    if (a.statsDir.empty())
      a.statsDir = dirName(replaceAll(a.base, "{N}", std::to_string(a.Ns[0]))) + "/../stats";

    std::string time = a.time;
    if (time == "latest") time = latestTime(a.caseDir + "/processor0");

    ConvertOptions copt;
    copt.field = a.field;
    copt.time = time;
    copt.qualityMode = a.qualityMode;
    copt.splitBadPyramids = a.splitBadPyramids;
    copt.seamDir = a.seamDir;
    copt.wallMaskLayers = a.wallMask;
    copt.wallMaskValue = a.wallMaskValue;
    if (a.seamPrepass && a.seamDir.empty()) usage("--seam-prepass requires --seam-dir");
    if (!a.seamDir.empty()) mkdirs(a.seamDir);

    const int rank   = getenv("SLURM_PROCID") ? atoi(getenv("SLURM_PROCID")) : 0;
    const int ntasks = getenv("SLURM_NTASKS") ? atoi(getenv("SLURM_NTASKS")) : 1;
    const int bBegin = a.blocksBegin >= 0 ? a.blocksBegin : 0;
    const int bEnd   = a.blocksEnd   >= 0 ? a.blocksEnd   : numBlocks;

    char host[256] = {0};
    gethostname(host, sizeof(host) - 1);
    printf("#f2u: rank %i/%i on %s: case %s, %i processors, field '%s', time '%s', "
           "grouping %s, block size %i, blocks [%i,%i), threads %i\n",
           rank, ntasks, host, a.caseDir.c_str(), P, a.field.c_str(), time.c_str(),
           a.grouping.c_str(), B, bBegin, bEnd, omp_get_max_threads());
    fflush(stdout);

    mkdirs(a.statsDir);
    for (int N : a.Ns) mkdirs(dirName(replaceAll(a.base, "{N}", std::to_string(N))));

    if (a.seamPrepass) {
      for (int b = bBegin + rank; b < bEnd; b += ntasks) {
        const double t0 = now();
        const int p0 = b * B;
        std::string firstError;
        int64_t pts = 0;
#pragma omp parallel for schedule(dynamic,1) reduction(+:pts)
        for (int i = 0; i < B; i++) {
          try {
            ProcStats st;
            writeSeamFile(a.caseDir, order[p0 + i], copt, st);
            pts += st.seamPoints;
          } catch (std::exception &e) {
#pragma omp critical
            if (firstError.empty()) firstError = e.what();
          }
        }
        if (!firstError.empty()) throw std::runtime_error(firstError);
        printf("#f2u: seam pre-pass block %i done: %lli processor-patch points (%.1f s)\n",
               b, (long long)pts, now() - t0);
        fflush(stdout);
      }
      return 0;
    }

    for (int b = bBegin + rank; b < bEnd; b += ntasks) {
      const double t0 = now();
      const int p0 = b * B, p1 = p0 + B;
      std::vector<ProcMesh> procs(B);
      std::string firstError;
      int done = 0;
#pragma omp parallel for schedule(dynamic,1)
      for (int i = 0; i < B; i++) {
        try {
          convertProcessor(a.caseDir, order[p0 + i], copt, procs[i]);
        } catch (std::exception &e) {
#pragma omp critical
          if (firstError.empty()) firstError = e.what();
        }
#pragma omp critical
        {
          done++;
          if (done % 8 == 0 || done == B) {
            printf("#f2u: block %i: converted %i/%i processors (%.1f s)\n", b, done, B, now() - t0);
            fflush(stdout);
          }
        }
      }
      if (!firstError.empty()) throw std::runtime_error(firstError);
      const double t1 = now();

      // ---- write parts for all requested N
      const std::string attrName = a.field;
      std::vector<PartRecord> records;
      for (int N : a.Ns) {
        const int ppp = P / N; // processors per part
        const int iBegin = p0 / ppp, iEnd = p1 / ppp;
        std::vector<PartRecord> recs(iEnd - iBegin);
#pragma omp parallel for schedule(dynamic,1)
        for (int i = iBegin; i < iEnd; i++) {
          try {
            std::vector<const ProcMesh *> list;
            for (int p = i * ppp; p < (i + 1) * ppp; p++) list.push_back(&procs[p - p0]);
            PartRecord r = writePart(list, partFileName(a.base, N, i), attrName);
            r.N = N; r.index = i;
            for (int p = i * ppp; p < (i + 1) * ppp; p++) r.procs.push_back(order[p]);
            recs[i - iBegin] = r;
          } catch (std::exception &e) {
#pragma omp critical
            if (firstError.empty()) firstError = e.what();
          }
        }
        if (!firstError.empty()) throw std::runtime_error(firstError);
        records.insert(records.end(), recs.begin(), recs.end());
      }
      const double t2 = now();

      // ---- block statistics (text, merged by --summarize)
      const std::string sname = a.statsDir + "/block_" + std::to_string(b);
      {
        std::ofstream o(sname + ".txt.tmp");
        o.precision(17);
        o << "block " << b << ' ' << p0 << ' ' << p1 << ' ' << P << ' ' << numBlocks << '\n';
        o << "field " << a.field << '\n';
        o << "time " << time << '\n';
        o << "timing " << (t1 - t0) << ' ' << (t2 - t1) << '\n';
        for (const auto &pm : procs) {
          o << "proc " << pm.procID << ' ';
          writeStats(o, pm.stats);
          o << ' ';
          writeBox(o, pm.bounds);
          o << '\n';
          for (const auto &pb : pm.patchBounds) {
            o << "patch " << pm.procID << ' ' << pb.first << ' '
              << pm.patchTypes.at(pb.first) << ' ' << pb.second.first << ' ';
            writeBox(o, pb.second.second);
            o << '\n';
          }
        }
        for (const auto &r : records) {
          o << "part " << r.N << ' ' << r.index << ' '
            << r.cells << ' ' << r.vertices << ' ' << r.tets << ' ' << r.pyrs << ' '
            << r.wedges << ' ' << r.hexes << ' ' << r.indices << ' ' << r.fileBytes << ' '
            << r.vmin << ' ' << r.vmax << ' ';
          writeBox(o, r.bounds);
          o << ' ' << r.file << ' ';
          for (size_t k = 0; k < r.procs.size(); k++) o << (k ? "," : "") << r.procs[k];
          o << '\n';
        }
      }
      // vertex samples (position + value) for percentiles, also per region
      {
        std::vector<float> samples;
        for (const auto &pm : procs)
          for (size_t i = 0; i < pm.values.size(); i += a.sampleStride) {
            const auto &v = pm.vertices[i];
            samples.insert(samples.end(), { v.x, v.y, v.z, pm.values[i] });
          }
        FILE *fp = fopen((sname + ".samples4.tmp").c_str(), "wb");
        if (!fp) throw std::runtime_error("cannot write samples");
        fwrite(samples.data(), sizeof(float), samples.size(), fp);
        fclose(fp);
      }
      rename((sname + ".samples4.tmp").c_str(), (sname + ".samples4").c_str());
      rename((sname + ".txt.tmp").c_str(), (sname + ".txt").c_str());

      ProcStats tot;
      size_t mem = 0;
      for (const auto &pm : procs) { tot.add(pm.stats); mem += pm.memoryBytes(); }
      printf("#f2u: block %i done: %lli cells -> %lli verts, %lli tets %lli pyrs %lli wedges "
             "%lli hexes; decomposed %lli poly + %lli quality + %lli topology; "
             "inverted t/p/w/h %lli/%lli/%lli/%lli; not-star t/p/w/h %lli/%lli/%lli/%lli; "
             "bad-corner hex/wedge %lli/%lli; quads split %lli; seam points %lli matched %lli; "
             "convert %.1f s, write %.1f s, mem %.2f GB\n",
             b, (long long)tot.cells, (long long)tot.outVertices, (long long)tot.outTets,
             (long long)tot.outPyrs, (long long)tot.outWedges, (long long)tot.outHexes,
             (long long)tot.decomposedPoly, (long long)tot.decomposedQuality,
             (long long)tot.decomposedTopology,
             (long long)tot.invertedTets, (long long)tot.invertedPyrs,
             (long long)tot.invertedWedges, (long long)tot.invertedHexes,
             (long long)tot.notStarTets, (long long)tot.notStarPyrs,
             (long long)tot.notStarWedges, (long long)tot.notStarHexes,
             (long long)tot.badCornerHexes, (long long)tot.badCornerWedges,
             (long long)tot.quadsSplitToTets, (long long)tot.seamPoints, (long long)tot.seamMatched,
             t1 - t0, t2 - t1, mem / 1e9);
      fflush(stdout);
    }
    return 0;
  }

  // ==================================================================
  // summarize mode
  // ==================================================================
  int runSummarize(Args &a)
  {
    if (a.statsDir.empty()) usage("--stats-dir missing");
    if (a.summary.empty()) usage("--summary missing");
    std::vector<std::string> files;
    {
      DIR *d = opendir(a.statsDir.c_str());
      if (!d) throw std::runtime_error("cannot open " + a.statsDir);
      while (dirent *e = readdir(d)) {
        std::string n = e->d_name;
        if (n.rfind("block_", 0) == 0 && n.size() > 4 && n.substr(n.size() - 4) == ".txt")
          files.push_back(a.statsDir + "/" + n);
      }
      closedir(d);
    }
    if (files.empty()) throw std::runtime_error("no block statistics in " + a.statsDir);

    ProcStats total;
    Box globalBounds;
    std::map<int, ProcStats> perProc;
    std::map<int, Box> procBounds;
    std::map<std::string, std::pair<int64_t,Box>> patches;
    std::map<std::string, std::string> patchTypes;
    std::map<int, std::map<int, PartRecord>> parts; // N -> index -> rec
    std::string field, time;
    int P = 0, numBlocks = 0;
    double maxConvert = 0, maxWrite = 0;
    std::vector<std::array<float,3>> samplePos;
    std::set<int> blocksSeen;
    std::vector<float> samples;

    for (const auto &fn : files) {
      std::ifstream in(fn);
      std::string line;
      while (std::getline(in, line)) {
        std::istringstream ls(line);
        std::string tag;
        ls >> tag;
        if (tag == "block") {
          int b, p0, p1; ls >> b >> p0 >> p1 >> P >> numBlocks; blocksSeen.insert(b);
        } else if (tag == "field") ls >> field;
        else if (tag == "time") ls >> time;
        else if (tag == "timing") {
          double c, w; ls >> c >> w;
          maxConvert = std::max(maxConvert, c); maxWrite = std::max(maxWrite, w);
        } else if (tag == "proc") {
          int id; ls >> id;
          ProcStats s; readStats(ls, s);
          Box bb; readBox(ls, bb);
          perProc[id] = s; procBounds[id] = bb;
          total.add(s);
          globalBounds.extend(bb);
        } else if (tag == "patch") {
          int id; std::string name, type; int64_t nf; Box bb;
          ls >> id >> name >> type >> nf; readBox(ls, bb);
          patches[name].first += nf;
          patches[name].second.extend(bb);
          patchTypes[name] = type;
        } else if (tag == "part") {
          PartRecord r;
          ls >> r.N >> r.index >> r.cells >> r.vertices
             >> r.tets >> r.pyrs >> r.wedges >> r.hexes >> r.indices >> r.fileBytes
             >> r.vmin >> r.vmax;
          readBox(ls, r.bounds);
          std::string plist;
          ls >> r.file >> plist;
          for (int pid : parseIntList(plist)) r.procs.push_back(pid);
          parts[r.N][r.index] = r;
        }
      }
      // samples: ".samples4" = (x,y,z,value) float quadruples; old
      // ".samples" files hold values only (no region statistics)
      const std::string base_ = fn.substr(0, fn.size() - 4);
      const int64_t sz4 = fileSize(base_ + ".samples4");
      const int64_t sz1 = fileSize(base_ + ".samples");
      if (sz4 > 0) {
        std::vector<float> buf(sz4 / sizeof(float));
        FILE *fp = fopen((base_ + ".samples4").c_str(), "rb");
        if (fp) { size_t r = fread(buf.data(), sizeof(float), buf.size(), fp); (void)r; fclose(fp); }
        for (size_t i = 0; i + 3 < buf.size(); i += 4) {
          samples.push_back(buf[i+3]);
          samplePos.push_back({ buf[i], buf[i+1], buf[i+2] });
        }
      } else if (sz1 > 0) {
        const size_t off = samples.size();
        samples.resize(off + sz1 / sizeof(float));
        FILE *fp = fopen((base_ + ".samples").c_str(), "rb");
        if (fp) { size_t r = fread(samples.data() + off, sizeof(float), sz1 / sizeof(float), fp); (void)r; fclose(fp); }
      }
    }
    if (!samplePos.empty() && samplePos.size() != samples.size())
      throw std::runtime_error("mixed old/new sample files in " + a.statsDir);

    // ---- value distributions (all samples + optional regions)
    struct Region { std::string name; double b[6]; };
    std::vector<Region> regions;
    for (const auto &rs : a.regions) {
      Region r;
      const size_t c = rs.find(':');
      if (c == std::string::npos) throw std::runtime_error("--region name:x0,x1,y0,y1,z0,z1");
      r.name = rs.substr(0, c);
      std::stringstream ss(rs.substr(c + 1));
      std::string t; int k = 0;
      while (std::getline(ss, t, ',') && k < 6) r.b[k++] = std::stod(t);
      if (k != 6) throw std::runtime_error("--region needs 6 numbers: " + rs);
      regions.push_back(r);
    }
    const std::vector<double> qs = { 0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 0.9999 };
    const std::vector<std::string> qn = { "p0.1", "p1", "p5", "p25", "p50", "p75", "p90", "p95", "p99", "p99.9", "p99.99" };
    std::ostringstream distJs, distTxt;
    auto describe = [&](const std::string &name, std::vector<float> v, const std::string &what) {
      std::sort(v.begin(), v.end());
      auto q = [&](const std::vector<float> &x, double f) -> double {
        if (x.empty()) return NAN;
        return x[std::min(x.size() - 1, (size_t)std::llround(f * (x.size() - 1)))];
      };
      std::vector<float> pos;
      for (float x : v) if (x > 0.f) pos.push_back(x);
      char buf[256];
      distJs << "    " << jstr(name) << ": {\"region\": " << jstr(what) << ", \"samples\": " << v.size()
             << ", \"min\": " << jnum(v.empty() ? NAN : v.front())
             << ", \"max\": " << jnum(v.empty() ? NAN : v.back())
             << ", \"fractionPositive\": " << jnum(v.empty() ? NAN : double(pos.size()) / v.size());
      snprintf(buf, sizeof(buf), "%-10s %-34s samples %zu, min %g, max %g, fraction > 0: %.4f\n",
               name.c_str(), what.c_str(), v.size(), v.empty() ? NAN : v.front(),
               v.empty() ? NAN : v.back(), v.empty() ? NAN : double(pos.size()) / v.size());
      distTxt << buf << "           all      :";
      distJs << ", \"percentiles\": {";
      for (size_t i = 0; i < qs.size(); i++) {
        distJs << (i ? ", " : "") << jstr(qn[i]) << ": " << jnum(q(v, qs[i]));
        snprintf(buf, sizeof(buf), " %s %.4g", qn[i].c_str(), q(v, qs[i]));
        distTxt << buf;
      }
      distJs << "}, \"percentilesOfPositive\": {";
      distTxt << "\n           values>0 :";
      for (size_t i = 0; i < qs.size(); i++) {
        distJs << (i ? ", " : "") << jstr(qn[i]) << ": " << jnum(q(pos, qs[i]));
        snprintf(buf, sizeof(buf), " %s %.4g", qn[i].c_str(), q(pos, qs[i]));
        distTxt << buf;
      }
      distJs << "}}";
      distTxt << "\n";
    };
    describe("all", samples, "whole domain");
    for (const auto &r : regions) {
      std::vector<float> v;
      for (size_t i = 0; i < samplePos.size(); i++) {
        const auto &P_ = samplePos[i];
        if (P_[0] >= r.b[0] && P_[0] <= r.b[1] && P_[1] >= r.b[2] && P_[1] <= r.b[3]
            && P_[2] >= r.b[4] && P_[2] <= r.b[5]) v.push_back(samples[i]);
      }
      char what[160];
      snprintf(what, sizeof(what), "x[%g,%g] y[%g,%g] z[%g,%g]",
               r.b[0], r.b[1], r.b[2], r.b[3], r.b[4], r.b[5]);
      distJs << ",\n";
      describe(r.name, v, what);
    }

    std::vector<std::string> problems;
    if ((int)blocksSeen.size() != numBlocks)
      problems.push_back("only " + std::to_string(blocksSeen.size()) + " of "
                         + std::to_string(numBlocks) + " blocks present");
    if ((int)perProc.size() != P)
      problems.push_back("only " + std::to_string(perProc.size()) + " of "
                         + std::to_string(P) + " processors present");

    // percentiles of vertex values (sampled)
    std::sort(samples.begin(), samples.end());
    auto pct = [&](double q) -> double {
      if (samples.empty()) return NAN;
      size_t i = (size_t)std::llround(q * (samples.size() - 1));
      return samples[std::min(i, samples.size() - 1)];
    };

    // value range over all parts of the first N (identical for every N)
    double vmin = +1e300, vmax = -1e300;
    for (auto &pn : parts) for (auto &pr : pn.second) {
      vmin = std::min(vmin, pr.second.vmin); vmax = std::max(vmax, pr.second.vmax);
    }

    // aircraft bounds: all wall patches, plus all non-far-field patches
    Box wallBounds, aircraftBounds;
    const std::set<std::string> farField = { "front", "back", "left", "right", "inlet", "outlet" };
    for (auto &p : patches) {
      if (patchTypes[p.first] == "wall") wallBounds.extend(p.second.second);
      if (!farField.count(p.first)) aircraftBounds.extend(p.second.second);
    }

    std::ostringstream js;
    js << "{\n";
    js << "  \"field\": " << jstr(field) << ",\n";
    js << "  \"time\": " << jstr(time) << ",\n";
    js << "  \"numProcessors\": " << P << ",\n";
    js << "  \"bounds\": " << jbox(globalBounds) << ",\n";
    js << "  \"wallPatchBounds\": " << jbox(wallBounds) << ",\n";
    js << "  \"aircraftBounds\": " << jbox(aircraftBounds) << ",\n";
    js << "  \"valueRange\": [" << jnum(vmin) << ", " << jnum(vmax) << "],\n";
    js << "  \"cellValueRange\": [" << jnum(total.cellMin) << ", " << jnum(total.cellMax) << "],\n";
    js << "  \"vertexValuePercentiles\": {\"samples\": " << samples.size()
       << ", \"p0.1\": " << jnum(pct(0.001)) << ", \"p1\": " << jnum(pct(0.01))
       << ", \"p5\": " << jnum(pct(0.05)) << ", \"p50\": " << jnum(pct(0.5))
       << ", \"p95\": " << jnum(pct(0.95)) << ", \"p99\": " << jnum(pct(0.99))
       << ", \"p99.9\": " << jnum(pct(0.999)) << "},\n";
    js << "  \"distributions\": {\n" << distJs.str() << "\n  },\n";
    js << "  \"input\": {\"cells\": " << total.cells << ", \"points\": " << total.points
       << ", \"faces\": " << total.faces << ", \"internalFaces\": " << total.internalFaces
       << ", \"triFaces\": " << total.facesBySize[0] << ", \"quadFaces\": " << total.facesBySize[1]
       << ", \"polyFaces\": " << total.facesBySize[2] << ", \"nanCells\": " << total.nanCells << "},\n";
    js << "  \"cellClassification\": {\"hex\": " << total.nativeHex
       << ", \"wedge\": " << total.nativeWedge << ", \"pyr\": " << total.nativePyr
       << ", \"tet\": " << total.nativeTet
       << ", \"decomposedPolyhedra\": " << total.decomposedPoly
       << ", \"decomposedBadQuality\": " << total.decomposedQuality
       << ", \"decomposedBadTopology\": " << total.decomposedTopology << "},\n";
    js << "  \"output\": {\"vertices\": " << total.outVertices
       << ", \"addedCellCentres\": " << total.addedCellCentres
       << ", \"addedFaceCentres\": " << total.addedFaceCentres
       << ", \"tets\": " << total.outTets << ", \"pyrs\": " << total.outPyrs
       << ", \"wedges\": " << total.outWedges << ", \"hexes\": " << total.outHexes
       << ", \"elements\": " << (total.outTets + total.outPyrs + total.outWedges + total.outHexes)
       << "},\n";
    js << "  \"elementChecks\": {\"invertedTets\": " << total.invertedTets
       << ", \"invertedPyrs\": " << total.invertedPyrs << ", \"invertedWedges\": " << total.invertedWedges
       << ", \"invertedHexes\": " << total.invertedHexes
       << ", \"notStarTets\": " << total.notStarTets << ", \"notStarPyrs\": " << total.notStarPyrs
       << ", \"notStarWedges\": " << total.notStarWedges << ", \"notStarHexes\": " << total.notStarHexes
       << ", \"badCornerHexes\": " << total.badCornerHexes
       << ", \"badCornerWedges\": " << total.badCornerWedges
       << ", \"quadsSplitToTets\": " << total.quadsSplitToTets << "},\n";
    js << "  \"wallMask\": {\"maskedCells\": " << total.maskedCells
       << ", \"maskedVertices\": " << total.maskedVertices << "},\n";
    js << "  \"seams\": {\"applied\": " << (total.seamContributions > 0 ? "true" : "false")
       << ", \"processorPatchPoints\": " << total.seamPoints
       << ", \"matched\": " << total.seamMatched
       << ", \"contributions\": " << total.seamContributions << "},\n";
    js << "  \"timing\": {\"maxBlockConvertSeconds\": " << jnum(maxConvert)
       << ", \"maxBlockWriteSeconds\": " << jnum(maxWrite) << "},\n";
    js << "  \"patches\": {\n";
    {
      size_t k = 0;
      for (auto &p : patches) {
        js << "    " << jstr(p.first) << ": {\"type\": " << jstr(patchTypes[p.first])
           << ", \"faces\": " << p.second.first << ", \"bounds\": " << jbox(p.second.second) << "}"
           << (++k < patches.size() ? ",\n" : "\n");
      }
    }
    js << "  },\n";

    // GPU memory estimate per part (barney via haystack/ANARI):
    //   vertices: float3 position + float scalar  = 16 B / vertex
    //   elements: uint32 index per element vertex + uint8 type + uint32 begin
    //   BVH/acceleration: ~64 B (low) .. ~160 B (high) per element
    auto gpuLow  = [](const PartRecord &r) {
      const int64_t e = r.tets + r.pyrs + r.wedges + r.hexes;
      return 16.0*r.vertices + 4.0*r.indices + 5.0*e + 64.0*e; };
    auto gpuHigh = [](const PartRecord &r) {
      const int64_t e = r.tets + r.pyrs + r.wedges + r.hexes;
      return 16.0*r.vertices + 4.0*r.indices + 5.0*e + 160.0*e; };

    js << "  \"partitions\": {\n";
    std::ostringstream txt;
    size_t kn = 0;
    for (auto &pn : parts) {
      const int N = pn.first;
      int64_t sv = 0, se = 0, sc = 0, sb = 0, maxE = 0, maxV = 0;
      double maxLow = 0, maxHigh = 0, boxVol = 0;
      auto bvol = [](const Box &b) {
        return b.empty() ? 0.0 : (b.hi[0]-b.lo[0])*(b.hi[1]-b.lo[1])*(b.hi[2]-b.lo[2]); };
      for (auto &pr : pn.second) {
        const PartRecord &r = pr.second;
        boxVol += bvol(r.bounds);
        const int64_t e = r.tets + r.pyrs + r.wedges + r.hexes;
        sv += r.vertices; se += e; sc += r.cells; sb += r.fileBytes;
        maxE = std::max(maxE, e); maxV = std::max(maxV, r.vertices);
        maxLow = std::max(maxLow, gpuLow(r)); maxHigh = std::max(maxHigh, gpuHigh(r));
      }
      if ((int)pn.second.size() != N)
        problems.push_back("N=" + std::to_string(N) + ": only " + std::to_string(pn.second.size())
                           + " parts present");
      if (sc != total.cells)
        problems.push_back("N=" + std::to_string(N) + ": sum of part cells != total cells");
      js << "    " << jstr(std::to_string(N)) << ": {\"parts\": " << pn.second.size()
         << ", \"vertices\": " << sv << ", \"elements\": " << se << ", \"cells\": " << sc
         << ", \"fileBytes\": " << sb << ", \"maxPartElements\": " << maxE
         << ", \"maxPartVertices\": " << maxV
         << ", \"partBoxOverlap\": " << jnum(boxVol / std::max(1e-300, bvol(globalBounds)))
         << ", \"gpuBytesPerPartEstimate\": [" << jnum(maxLow) << ", " << jnum(maxHigh) << "],\n";
      js << "      \"list\": [\n";
      size_t k = 0;
      for (auto &pr : pn.second) {
        const PartRecord &r = pr.second;
        js << "        {\"index\": " << r.index << ", \"file\": " << jstr(r.file)
           << ", \"processors\": [" << [&]{ std::string t; for (size_t k = 0; k < r.procs.size(); k++)
                                                 t += (k ? ", " : "") + std::to_string(r.procs[k]);
                                               return t; }() << "]"
           << ", \"cells\": " << r.cells << ", \"vertices\": " << r.vertices
           << ", \"tets\": " << r.tets << ", \"pyrs\": " << r.pyrs
           << ", \"wedges\": " << r.wedges << ", \"hexes\": " << r.hexes
           << ", \"fileBytes\": " << r.fileBytes
           << ", \"valueRange\": [" << jnum(r.vmin) << ", " << jnum(r.vmax) << "]"
           << ", \"bounds\": " << jbox(r.bounds) << "}"
           << (++k < pn.second.size() ? ",\n" : "\n");
      }
      js << "      ]}" << (++kn < parts.size() ? ",\n" : "\n");

      char buf[512];
      snprintf(buf, sizeof(buf),
               "N=%-4d parts=%zu  vertices=%lli  elements=%lli  cells=%lli  files=%.1f GB  "
               "max/part: %lli elements, %lli vertices, est. GPU %.2f-%.2f GB; "
               "sum(part boxes)/domain %.2f\n",
               N, pn.second.size(), (long long)sv, (long long)se, (long long)sc, sb / 1e9,
               (long long)maxE, (long long)maxV, maxLow / 1e9, maxHigh / 1e9,
               boxVol / std::max(1e-300, bvol(globalBounds)));
      txt << buf;
    }
    js << "  },\n";
    js << "  \"problems\": [";
    for (size_t i = 0; i < problems.size(); i++) js << (i ? ", " : "") << jstr(problems[i]);
    js << "]\n}\n";

    {
      std::ofstream o(a.summary + ".json");
      o << js.str();
    }
    {
      std::ofstream o(a.summary + ".txt");
      char buf[1024];
      o << "foam2umesh summary: field " << field << ", time " << time << ", "
        << P << " processors\n";
      snprintf(buf, sizeof(buf),
               "bounds          : (%g %g %g) - (%g %g %g)\n"
               "wall bounds     : (%g %g %g) - (%g %g %g)\n"
               "aircraft bounds : (%g %g %g) - (%g %g %g)\n"
               "value range     : [%g, %g] (cells [%g, %g])\n"
               "percentiles     : p0.1 %g  p1 %g  p5 %g  p50 %g  p95 %g  p99 %g  p99.9 %g\n",
               globalBounds.lo[0], globalBounds.lo[1], globalBounds.lo[2],
               globalBounds.hi[0], globalBounds.hi[1], globalBounds.hi[2],
               wallBounds.lo[0], wallBounds.lo[1], wallBounds.lo[2],
               wallBounds.hi[0], wallBounds.hi[1], wallBounds.hi[2],
               aircraftBounds.lo[0], aircraftBounds.lo[1], aircraftBounds.lo[2],
               aircraftBounds.hi[0], aircraftBounds.hi[1], aircraftBounds.hi[2],
               vmin, vmax, total.cellMin, total.cellMax,
               pct(0.001), pct(0.01), pct(0.05), pct(0.5), pct(0.95), pct(0.99), pct(0.999));
      o << buf;
      snprintf(buf, sizeof(buf), "wall mask       : %lli cells, %lli vertices set to the mask value\n",
               (long long)total.maskedCells, (long long)total.maskedVertices);
      o << buf;
      o << "sampled vertex-value distributions (every k-th vertex):\n" << distTxt.str();
      snprintf(buf, sizeof(buf),
               "input           : %lli cells, %lli points, %lli faces (%lli tri, %lli quad, %lli >4)\n"
               "native cells    : %lli hex, %lli wedge, %lli pyr, %lli tet\n"
               "decomposed cells: %lli polyhedra, %lli bad-quality, %lli bad-topology\n"
               "output          : %lli vertices (+%lli cell centres, +%lli face centres)\n"
               "                  %lli tets, %lli pyrs, %lli wedges, %lli hexes = %lli elements\n"
               "inverted (vol<=0): %lli tets, %lli pyrs, %lli wedges, %lli hexes\n"
               "not star-shaped : %lli tets, %lli pyrs, %lli wedges, %lli hexes\n"
               "bad corner Jac. : %lli hexes, %lli wedges (native cells, before decision)\n"
               "quads -> 2 tets : %lli\n"
               "seams           : %lli processor-patch points, %lli matched on other processors (%lli contributions)\n"
               "NaN cells       : %lli\n"
               "timing          : max block convert %.1f s, write %.1f s\n",
               (long long)total.cells, (long long)total.points, (long long)total.faces,
               (long long)total.facesBySize[0], (long long)total.facesBySize[1],
               (long long)total.facesBySize[2],
               (long long)total.nativeHex, (long long)total.nativeWedge,
               (long long)total.nativePyr, (long long)total.nativeTet,
               (long long)total.decomposedPoly, (long long)total.decomposedQuality,
               (long long)total.decomposedTopology,
               (long long)total.outVertices, (long long)total.addedCellCentres,
               (long long)total.addedFaceCentres,
               (long long)total.outTets, (long long)total.outPyrs, (long long)total.outWedges,
               (long long)total.outHexes,
               (long long)(total.outTets + total.outPyrs + total.outWedges + total.outHexes),
               (long long)total.invertedTets, (long long)total.invertedPyrs,
               (long long)total.invertedWedges, (long long)total.invertedHexes,
               (long long)total.notStarTets, (long long)total.notStarPyrs,
               (long long)total.notStarWedges, (long long)total.notStarHexes,
               (long long)total.badCornerHexes, (long long)total.badCornerWedges,
               (long long)total.quadsSplitToTets,
               (long long)total.seamPoints, (long long)total.seamMatched,
               (long long)total.seamContributions,
               (long long)total.nanCells, maxConvert, maxWrite);
      o << buf;
      o << txt.str();
      o << "patches (non-empty, non-processor):\n";
      for (auto &p : patches) {
        const Box &b = p.second.second;
        snprintf(buf, sizeof(buf), "  %-28s %-14s %10lli faces  (%g %g %g) - (%g %g %g)\n",
                 p.first.c_str(), patchTypes[p.first].c_str(), (long long)p.second.first,
                 b.lo[0], b.lo[1], b.lo[2], b.hi[0], b.hi[1], b.hi[2]);
        o << buf;
      }
      for (auto &pr : problems) o << "PROBLEM: " << pr << "\n";
    }
    std::ifstream t(a.summary + ".txt");
    std::cout << t.rdbuf();
    return problems.empty() ? 0 : 2;
  }

  // ==================================================================
  // check mode
  // ==================================================================
  int runCheck(const std::string &fileName)
  {
    const double t0 = now();
    umesh::UMesh::SP mesh = umesh::UMesh::loadFrom(fileName);
    const double t1 = now();
    std::cout << "#f2u.check: loaded " << fileName << " in " << (t1 - t0) << " s\n";
    std::cout << "#f2u.check: " << mesh->toString(false) << "\n";
    const umesh::box3f b = mesh->getBounds();
    const umesh::range1f r = mesh->getValueRange();
    printf("#f2u.check: bounds (%g %g %g) - (%g %g %g), value range [%g, %g], attribute '%s'\n",
           b.lower.x, b.lower.y, b.lower.z, b.upper.x, b.upper.y, b.upper.z,
           r.lower, r.upper, mesh->perVertex ? mesh->perVertex->name.c_str() : "<none>");
    int64_t errors = 0;
    const int64_t nv = mesh->vertices.size();
    if (!mesh->perVertex || (int64_t)mesh->perVertex->values.size() != nv) {
      printf("#f2u.check: ERROR per-vertex attribute size mismatch\n"); errors++;
    }
    int64_t nanV = 0, nanS = 0;
    for (auto &v : mesh->vertices)
      if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z)) nanV++;
    if (mesh->perVertex) for (float s : mesh->perVertex->values) if (!std::isfinite(s)) nanS++;
    printf("#f2u.check: non-finite vertices %lli, non-finite scalars %lli\n",
           (long long)nanV, (long long)nanS);
    errors += nanV + nanS;

    std::vector<char> used(nv, 0);
    int64_t badIdx = 0;
    auto chk = [&](const int *idx, int n) {
      for (int k = 0; k < n; k++) {
        if (idx[k] < 0 || idx[k] >= nv) { badIdx++; return false; }
        used[idx[k]] = 1;
      }
      return true;
    };
    // (a) orientation / star test of this converter (checkElement)
    // (b) umesh's own orientation tests from apps/fixNegativeVolumeElements.cpp
    //     (elements that umesh would "fix" by swapping)
    // (c) barney's tet test (tetScalar: all four implicit-plane distances
    //     of the centroid must be >= 0)
    int64_t inv[4] = {0,0,0,0}, star[4] = {0,0,0,0}, umeshSwap[4] = {0,0,0,0}, barneyTet = 0;
    double minVol[4] = {1e300,1e300,1e300,1e300};
    umesh::vec3f v[8];
    auto avg = [&](std::initializer_list<int> ids) {
      double c[3] = {0,0,0};
      for (int i : ids) for (int k = 0; k < 3; k++) c[k] += v[i][k];
      std::array<double,3> r = { c[0]/ids.size(), c[1]/ids.size(), c[2]/ids.size() };
      return r;
    };
    auto account = [&](int type, const VolCheck &vc) {
      if (!(vc.total > 0.0)) inv[type]++;
      if (!vc.allPositive) star[type]++;
      minVol[type] = std::min(minVol[type], vc.total);
    };
    for (auto &e : mesh->tets) {
      if (!chk(&e.x, 4)) continue;
      for (int k = 0; k < 4; k++) v[k] = mesh->vertices[e[k]];
      account(0, checkElement(v, tetFaces));
      if (tetVol(v[0],v[1],v[2],v[3]) < 0.f) umeshSwap[0]++;
      // barney: evalToImplicitPlane(P,a,b,c) = dot(P-a, cross(b-a,c-a)) at the centroid
      const auto P = avg({0,1,2,3});
      if (tetVol(v[0],v[1],v[2],P) < 0 || tetVol(v[0],v[3],v[1],P) < 0 ||
          tetVol(v[0],v[2],v[3],P) < 0 || tetVol(v[1],v[3],v[2],P) < 0) barneyTet++;
    }
    for (auto &e : mesh->pyrs) {
      if (!chk(&e.base.x, 5)) continue;
      for (int k = 0; k < 5; k++) v[k] = mesh->vertices[e[k]];
      account(1, checkElement(v, pyrFaces));
      if (tetVol(v[0],v[1],avg({0,1,2,3}),v[4]) < 0) umeshSwap[1]++;
    }
    for (auto &e : mesh->wedges) {
      if (!chk(&e.front.x, 6)) continue;
      for (int k = 0; k < 6; k++) v[k] = mesh->vertices[e[k]];
      account(2, checkElement(v, wedgeFaces));
      if (tetVol(v[3],v[4],v[5],avg({0,1,3,4})) < 0) umeshSwap[2]++;
    }
    for (auto &e : mesh->hexes) {
      if (!chk(&e.base.x, 8)) continue;
      for (int k = 0; k < 8; k++) v[k] = mesh->vertices[e[k]];
      account(3, checkElement(v, hexFaces));
      if (tetVol(v[0],v[1],avg({0,1,2,3}),avg({0,1,2,3,4,5,6,7})) < 0) umeshSwap[3]++;
    }
    int64_t unused = 0;
    for (char u : used) if (!u) unused++;

    // seams: vertices with bitwise identical coordinates (duplicated
    // points on processor boundaries inside this part) should carry
    // identical values when the seam fix was applied
    if (mesh->perVertex) {
      std::vector<int64_t> order(nv);
      for (int64_t i = 0; i < nv; i++) order[i] = i;
      const auto &V = mesh->vertices;
      std::sort(order.begin(), order.end(), [&](int64_t a, int64_t b) {
        if (V[a].x != V[b].x) return V[a].x < V[b].x;
        if (V[a].y != V[b].y) return V[a].y < V[b].y;
        return V[a].z < V[b].z; });
      const auto &S = mesh->perVertex->values;
      int64_t dupGroups = 0, dupDiffer = 0;
      double maxDiff = 0.0, sumRel = 0.0;
      const double range = std::max(1e-30, double(r.upper) - double(r.lower));
      for (int64_t i = 0; i < nv; ) {
        int64_t j = i + 1;
        while (j < nv && V[order[j]].x == V[order[i]].x && V[order[j]].y == V[order[i]].y
               && V[order[j]].z == V[order[i]].z) j++;
        if (j - i > 1) {
          dupGroups++;
          float lo = S[order[i]], hi = lo;
          for (int64_t k = i; k < j; k++) { lo = std::min(lo, S[order[k]]); hi = std::max(hi, S[order[k]]); }
          if (hi > lo) { dupDiffer++; maxDiff = std::max(maxDiff, double(hi - lo)); sumRel += (hi - lo) / range; }
        }
        i = j;
      }
      printf("#f2u.check: seams: %lli duplicated vertex positions, %lli with differing values "
             "(max diff %g, mean diff/range %g)\n", (long long)dupGroups, (long long)dupDiffer,
             maxDiff, dupDiffer ? sumRel / dupDiffer : 0.0);
    }
    printf("#f2u.check: out-of-range indices %lli, unreferenced vertices %lli\n",
           (long long)badIdx, (long long)unused);
    const char *names[4] = { "tets", "pyrs", "wedges", "hexes" };
    const size_t counts[4] = { mesh->tets.size(), mesh->pyrs.size(), mesh->wedges.size(), mesh->hexes.size() };
    for (int t = 0; t < 4; t++)
      printf("#f2u.check: %-6s %10zu: inverted (signed vol <= 0) %lli, not star-shaped %lli, "
             "umesh-fixNegativeVolume would swap %lli, min signed vol %g\n",
             names[t], counts[t], (long long)inv[t], (long long)star[t], (long long)umeshSwap[t],
             counts[t] ? minVol[t] : 0.0);
    printf("#f2u.check: tets failing barney's tetScalar() inside test at their centroid: %lli\n",
           (long long)barneyTet);
    errors += inv[0] + inv[1] + inv[2] + inv[3];
    errors += badIdx;
    printf("#f2u.check: %s (%.1f s)\n", errors ? "FAILED" : "OK", now() - t0);
    return errors ? 1 : 0;
  }

} // anonymous

int main(int ac, char **av)
{
  Args a;
  try {
    for (int i = 1; i < ac; i++) {
      const std::string arg = av[i];
      auto next = [&]() -> std::string {
        if (i + 1 >= ac) usage("missing value for " + arg);
        return av[++i];
      };
      if (arg == "-h" || arg == "--help") usage();
      else if (arg == "--case") a.caseDir = next();
      else if (arg == "-n") a.Ns = parseIntList(next());
      else if (arg == "-o") a.base = next();
      else if (arg == "--field") a.field = next();
      else if (arg == "--time") a.time = next();
      else if (arg == "--nprocs") a.nprocs = std::stoi(next());
      else if (arg == "--threads") a.threads = std::stoi(next());
      else if (arg == "--stats-dir") a.statsDir = next();
      else if (arg == "--blocks-begin") a.blocksBegin = std::stoi(next());
      else if (arg == "--blocks-end") a.blocksEnd = std::stoi(next());
      else if (arg == "--sample-stride") a.sampleStride = std::max(1, std::stoi(next()));
      else if (arg == "--no-quality-decompose") a.qualityMode = 0;
      else if (arg == "--quality") a.qualityMode = std::stoi(next());
      else if (arg == "--no-split-pyramids") a.splitBadPyramids = false;
      else if (arg == "--summarize") a.summarize = true;
      else if (arg == "--seam-dir") a.seamDir = next();
      else if (arg == "--seam-prepass") a.seamPrepass = true;
      else if (arg == "--grouping") a.grouping = next();
      else if (arg == "--region") a.regions.push_back(next());
      else if (arg == "--wall-mask") a.wallMask = std::stoi(next());
      else if (arg == "--wall-mask-value") a.wallMaskValue = std::stof(next());
      else if (arg == "--summary") a.summary = next();
      else if (arg == "--check") a.check = next();
      else usage("unknown argument '" + arg + "'");
    }
    if (!a.check.empty()) return runCheck(a.check);
    if (a.summarize) return runSummarize(a);
    return runConvert(a);
  } catch (std::exception &e) {
    fprintf(stderr, "#f2u: fatal error: %s\n", e.what());
    return 1;
  }
}
