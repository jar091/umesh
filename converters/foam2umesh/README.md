# foam2umesh

Converts a **decomposed OpenFOAM case** (binary `processorN/constant/polyMesh`
plus cell-centred `volScalarField`/`volVectorField` of one time step) into
**N stand-alone `.umesh` part files** for data-parallel rendering with
haystack/barney (`hsOffline part_0000.umesh ... part_<N-1>.umesh -ndg N`).

Written for the AIRBUS4 case (OpenFOAM v2506, snappyHexMesh, 1024 processors,
1,025,949,906 cells) but is generic for binary `label=32` cases.

## Build

Built as part of the umesh project (`converters/CMakeLists.txt`, C++17, requires
OpenMP; skipped with a message if OpenMP is not found). The `foam2umesh`
executable ends up in the build directory. With
`-DFOAM2UMESH_STATIC_RUNTIME=ON` libstdc++/libgcc are linked statically, so the
binary runs without loading the GCC module.

## Usage

```bash
foam2umesh --case <case> -n 128,64 -o <out>/n{N}/airbus4_p_t1 \
           [--field p] [--time latest|<dir>] [--grouping hierarchical:8,8,16] \
           [--seam-dir <dir>] [--stats-dir <dir>] [--blocks-begin b --blocks-end e]
foam2umesh --case <case> -n ... --seam-prepass --seam-dir <dir>      # before the conversion
foam2umesh --summarize --stats-dir <dir> --summary <prefix>          # -> <prefix>.json/.txt
foam2umesh --check <part.umesh>                                      # reload + validate
```

* Part `i` of `N` is written to `<base>_<iiii>.umesh` (`{N}` in `-o` is replaced).
  All requested `N` are written from one conversion pass.
* `--field`: scalars are used directly, vectors as magnitude (`U` -> |U|).
  The umesh attribute is named after the field. Default `p`, time `latest`.
* Parallelism: the processors are cut into `min(N)` *blocks*; each Slurm task
  (`SLURM_PROCID`/`SLURM_NTASKS`, no MPI needed) converts its blocks with OpenMP
  (one processor per thread) and writes all parts inside its blocks.
  A full run is: seam pre-pass, conversion, summary, validation (`--check`);
  for AIRBUS4 on Karolina CPU nodes ~2-3 min in total.

## Cell reconstruction

For every cell the faces are collected with outward orientation (OpenFOAM face
normals point out of `owner` and into `neighbour`, so neighbour faces are
reversed). Cells are classified by their face set:

| faces                 | element | vertex order written (umesh = VTK = barney)                         |
|-----------------------|---------|----------------------------------------------------------------------|
| 6 quads, 8 vertices   | HEX     | base = a face reversed (normal into the cell), top = partners along the 4 lateral edges |
| 2 tris + 3 quads      | WEDGE   | (v0,v1,v2) = a triangle in *outward* orientation, (v3,v4,v5) = lateral partners |
| 1 quad + 4 tris       | PYR     | base = quad reversed (normal towards apex), apex                     |
| 4 tris                | TET     | outward (a,b,c) + d -> (a,c,b,d)                                     |

The topology is verified (unique vertex count, every lateral face has exactly
two base vertices, partners form the opposite face); anything else is a
general polyhedron and is **decomposed**: a cell-centre vertex is added
(volume-weighted centroid computed like OpenFOAM's `primitiveMesh`), every quad
face becomes a pyramid and every triangle a tet towards the centre, and faces
with > 4 vertices get a (shared, per face) face-centre vertex and are fanned
into triangles -> tets. Shared face-centre vertices keep the mesh conforming
between neighbouring decomposed cells.

Quality (`--quality`, default 1): native hex/wedge/pyr cells that are not
star-shaped w.r.t. their centre (a folded face) are decomposed as well; in
decomposed cells a pyramid on a warped quad that is not star-shaped is replaced
by two positive tets when possible. Mode 2 additionally decomposes hexes/wedges
with a non-positive corner Jacobian, mode 0 never decomposes native shapes.

### Vertex-order convention (verified)

* **Tet** `(v0,v1,v2,v3)`: `dot(v3-v0, cross(v1-v0, v2-v0)) > 0`. barney's
  `tetScalar()` (`barney/native/umesh/common/UMeshField.h`) only accepts points
  with all four `evalToImplicitPlane()` values >= 0, so inverted tets would be
  invisible. Same convention as umesh `tetrahedralize.cpp` / `fixNegativeVolumeElements.cpp`.
* **Pyr** `(v0..v3, v4)`: base normal (right-hand rule) points to the apex.
* **Wedge** `(v0,v1,v2, v3,v4,v5)`: normal of `(v0,v1,v2)` points *away* from
  `(v3,v4,v5)` (VTK_WEDGE; umesh `FaceConn.cpp` lists `(v0,v2,v1)` as inward facet,
  `fixNegativeVolumeElements.cpp` swaps wedges with `vol(v3,v4,v5,centre(0,1,4,3)) < 0`).
* **Hex** `(v0..v3, v4..v7)`: base normal points to the top, `vi-v(i+4)` edges.

barney evaluates pyramids/wedges/hexes with the VTK/OpenVKL parametric Newton
iteration (`ElementIntersection.h`), which needs this vertex *topology* (its
sign is irrelevant); haystack passes the umesh arrays unchanged as VTK cell types
10/14/13/12.

Verification (`--check` and the converter statistics): every element is split
into face-fan sub-tets towards its centre (inward faces from umesh
`FaceConn.cpp`); signed volume > 0 proves the orientation, all sub-tets > 0 means
star-shaped. `--check` additionally runs umesh's own `fixNegativeVolumeElements`
criteria and barney's tet inside-test at the tet centroid. AIRBUS4: 0 inverted,
0 non-star elements, 0 elements umesh would swap, 0 tets failing barney's test.

## Scalars

Cell -> vertex: every original point gets the average of all incident cells.
Added cell-centre vertices get the cell value, face-centre vertices the average
of the face's vertex values.

**Seams**: with `--seam-dir` (recommended) the pre-pass writes, for
every processor, the points on its processor patches with (sum, count) of the
incident cell values; the conversion adds the contributions of all other
processors in the 2-ring of processor neighbours that have the bitwise-identical
point. Duplicated boundary vertices of neighbouring processors/parts therefore
carry identical values (AIRBUS4: all processor-patch points matched). Without
the seam pass vertex values on processor borders only average the cells of one
processor (small value seams, geometry is always crack-free because the
coordinates are identical).

## Partitioning

Parts are groups of `P/N` processors (every part a standalone umesh, vertices
local to the part). `--grouping contiguous` uses processor-id ranges.
AIRBUS4 was decomposed with `hierarchical (8 8 16)`; its processor id is
`xi + 8*(yi + 8*zi)`, so contiguous id ranges run across the x-slabs whose y/z
cuts do not line up (part bounding boxes overlapping 4-7x). Use
`--grouping hierarchical:8,8,16` (recommended): parts are consecutive
z-boxes of one (x,y) column (N=128: half column, N=64: full column, N=32/16:
2/4 neighbouring columns of one x-slab), giving nearly disjoint boxes
(sum of part boxes / domain ~1.1).

## Output / statistics

`<out>/<tag>_summary.json|txt`: global bounds, wall/aircraft patch bounds,
per-patch bounds, value range, vertex-value percentiles (every 128th vertex),
cell classification, element totals, check counters, per-N part list (files,
processors, counts, bounds, value ranges), GPU memory estimate.

## Limitations

* Binary format, 32-bit labels, `LSB` only; `internalField nonuniform` only
  (the `0/` directory with `uniform` fields cannot be converted).
* Vertex indices are 32 bit per part (fine up to ~2 G vertices per part).
* Points shared only with processors outside the 2-ring of face neighbours would
  miss those contributions (does not occur for this decomposition: 100% matched).
