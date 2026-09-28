# bhac2umesh

Converts a BHAC (MPI-AMRVAC based GRMHD code) snapshot `dataNNNN.dat` on a
spherical modified-Kerr-Schild (MKS) grid into unstructured meshes (`.umesh`)
that haystack/barney render, split into N spatially coherent parts for
data-parallel (MPI) rendering.

```
bhac2umesh <dataNNNN.dat> --par <amrvac.par> -o <outBase> [-n N] [--var rho] [--log|--lin]
           [--rmin r] [--rmax r] [--mks-h 0.35] [--mks-r0 0] [--spin a] [--mass M] [--gamma g]
bhacUmeshCheck <file.umesh> [...]
```

Output: `<outBase>_<i>.umesh` (i = 0..N-1) and `<outBase>.json` (bounds, value
range and percentiles, cell/vertex counts per part, snapshot time, ...).

Built as part of the umesh project (`converters/CMakeLists.txt`); the
executables `bhac2umesh`, `bhacUmeshCheck` and `barneyVolumeEmu` end up in the
build directory. OpenMP is optional (only speeds up the diagnostic modes).

## Input

* File layout, `.par` parsing and the forest walk are adapted from the
  SpaceConverter BHAC reader (`space-converter/src/bhac/bhac_extract_iolib.cpp`):
  `nleafs` blocks of `w(nx1,nx2,nx3,nw)` doubles (no ghost cells, all `nw`
  variables incl. the auxiliaries `lfac`, `xi`) + staggered `ws(0:nx1,0:nx2,0:nx3,nws)`,
  then the forest (one int32 per node), then the tail (`nx`, `eqpar`, `nleafs`,
  `levmax`, `ndim`, `ndir`, `nw`, `nws`, `neqpar`, `it`, `t`). Level-1 blocks
  are in ig1-fastest order (`SFC_ON_BASE` undefined in the torus build).
* The base grid (`nxlone`, `xprobmin/max`, `typeaxial`, `wnames`) comes from
  the `.par` file; the MKS parameters `h` and `R0` (coordpar in `amrvacusr.t`)
  are not stored anywhere in the output and are passed as options.
* Only single-level grids (`mxnest = 1`) are supported (AMR levels would need
  hanging nodes).

## Geometry

* Corner (x1, x2, x3) -> r = R0 + exp(x1), theta = x2 + h/2 sin(2 x2), phi = x3
  (BHAC `mod_metric.f`) -> x = r sin(theta) cos(phi), y = r sin(theta) sin(phi),
  z = r cos(theta) (flat embedding of the Kerr-Schild coordinates, units of M).
* Vertices are shared between cells, across the phi seam (phi = 0 == 2 pi)
  and on the axis (one vertex per radius and pole). Cells touching the axis
  become wedges, all others hexes.
* Vertex order = umesh = VTK (haystack passes the umesh arrays unchanged to
  ANARI/barney as VTK cell types 12/13):
  * hex: base (0,1,2,3) = (i,j),(i+1,j),(i+1,j+1),(i,j+1) at phi_k, top (4..7)
    the same at phi_k+1 -> base normal points to the top (positive Jacobian).
  * wedge: front triangle (0,1,2) on r_i, back (3,4,5) on r_i+1 with 3 over 0
    etc.; the right-hand normal of (0,1,2) points away from (3,4,5) - the
    criterion of `umesh/apps/fixNegativeVolumeElements.cpp`. (In the
    shape-function parametrisation that barney uses for prisms this is a
    negative Jacobian, which barney's Newton inversion does not care about.)
* The converter checks the umesh orientation criterion for every element and
  sums the signed volumes (ratio to the spherical shell volume ~0.999, the rest
  is the faceting of the spheres). `bhacUmeshCheck` reloads the files with
  `UMesh::loadFrom`, checks indices/NaNs/orientation/volumes and runs a copy of
  barney's point-in-hex/prism Newton test on points inside every spherical cell.

## Scalar

One per-vertex scalar (barney renders only per-vertex data): the average of
the (already transformed, e.g. log10) values of all cells sharing the vertex.

| `--var` | meaning |
|---|---|
| `rho` (default) | rest-mass density d / lfac (BHAC: D = Gamma rho, not densitized by sqrt(gamma)) |
| `p` | pressure (gamma-1)/gamma (xi / lfac^2 - rho), xi = Gamma^2 rho h |
| `temp` | p / rho |
| `bsq` | b^2 = B^2 / Gamma^2 + (S_i B^i / xi)^2, B^2 with the Kerr-Schild 3-metric in MKS coordinates |
| `sigma` | b^2 / rho |
| `beta` | 2 p / b^2 |
| any wname | stored variable as is (d s1 s2 s3 tau b1 b2 b3 Ds dtr1 lfac xi) |

Default scale: log10 for the derived variables, linear for stored ones
(override with `--log` / `--lin`). With log10, non-positive values get the
smallest positive value of the field.

## Partitioning

Recursive bisection of the logical cell box (after the optional radial cut)
along the direction with the largest angular extent (d ln r, d theta,
sin(theta) d phi), balanced by cell count. Every part is a spherical shell
sector with its own vertices (vertices on part borders are duplicated with
identical values). N need not be a power of two.

## Diagnostics and variants

* `--tets`: write tetrahedra instead of hexes/wedges (6 tets per hex around
  the 0-6 diagonal, 3 per axis cell; conforming between neighbours and across
  the phi seam; positive orientation as required by barney's tet test).
* `--analytic r|z|torus`: same grid, partitioning and cell->vertex averaging,
  but an analytic cell value (torus = log10 of a Gaussian torus around
  R_cyl = 20, sigma = 4). `bhacUmeshCheck --expect <f>` then compares the value
  that barney's shape functions interpolate at interior points with the exact
  field (`--bad-order-control` swaps one corner pair as a negative control).
* `--ref-render out.ppm [--ref-iso v --ref-camera ... --ref-fovy --ref-res W H]`:
  CPU ray-cast of the iso-surface of the per-vertex field (trilinear in the
  logical cell), no .umesh output.
* `barneyVolumeEmu -xf tf.xf -o prefix parts.umesh...`: CPU emulation (float) of
  barney's unstructured-volume path (MC grid + majorants from element boxes,
  DDA, Woodcock with fastLog, hex/prism Newton) for primary rays over all parts;
  `--global-majorant` replaces the MC grid by one majorant per part.
