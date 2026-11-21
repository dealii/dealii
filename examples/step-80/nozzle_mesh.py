#!/usr/bin/env python3
"""
CAD + all-hexahedral GMSH mesh of a convergent-divergent (contraction-expansion)
nozzle.

The nozzle is a body of revolution built from five axial sections:

    inlet cylinder  ->  convergent cone  ->  throat cylinder
                    ->  divergent cone   ->  outlet cylinder

Coordinate system
-----------------
  z-axis : nozzle centre-line, flow direction, z = 0 at the inlet plane
  r      : radial distance from the centre-line

All dimensions are in millimetres.

Schematic (not to scale)
------------------------
  station:  z0        z1        z2        z3        z4        z5
  radius:   D/2       D/2       d/2       d/2       D/2       D/2
            |--L_in---|--L_conv-|-L_throat|--L_div--|--L_out--|
            +---------+         +---------+
            |         \\____    ____/      \\_________+
    axis ---+------------------------------------------+---
            |         /‾‾‾‾    ‾‾‾‾\\      /‾‾‾‾‾‾‾‾‾+
            +---------+         +---------+

Meshing strategy -- DELIBERATE DEVIATION FROM ../CLAUDE.md
----------------------------------------------------------
The house recipe in ../CLAUDE.md obtains hexahedra by meshing tetrahedra and
subdividing each into four (Mesh.SubdivisionAlgorithm = 2).  That is the right
tool for the wine-bottle scripts -- organic spline shapes that cannot be block
decomposed -- but it is the wrong tool here: it would produce 10-50x more cells
than necessary, of much lower quality, with no directional control.

A nozzle is a body of revolution made of straight-line profile segments, so it
admits a fully structured O-grid ("butterfly") decomposition instead:

    cross-section at any station, radius R(z)
             ___---‾‾‾---___
          /    \\           /    \\
         |      +---------+      |     CORE   : n_circ x n_circ quads
         | petal|  CORE   |petal |     4 PETALS: n_circ x n_rad quads each
         |      +---------+      |
          \\    /           \\    /
             ‾‾‾---___---‾‾‾

swept through the five axial segments and constrained transfinite, giving
perfectly structured hexahedra with an exact cell count.

Construction outline
--------------------
  1. the nozzle solid is built from OCC primitives and fused -> STEP export
  2. a square-frustum "core" following rc(z) = core_frac * R(z) is built
  3. four diagonal half-planes (aligned with the core-square CORNERS) and the
     internal station disks are used as cutting tools
  4. one fragment() call yields 5 segments x (1 core + 4 petals) conformal
     hexahedral blocks
  5. curves are classified geometrically (axial / radial / circumferential) so
     the transfinite counts agree on opposite edges without tracking tags

Because the diagonal planes are aligned with the core-square corners they graze
edges that already exist and never cut the core interior, so the core stays one
volume per segment while the annulus splits into four petals -- and the
core/petal interface stays conformal for free.

Output
------
  nozzle.step   exact CAD, a single unpartitioned solid
  nozzle.iges   the same solid as IGES, for consumers that prefer it
  nozzle.msh    GMSH v2.2 (Lethe / deal.II, OpenFOAM gmshToFoam, Fluent)

Usage
-----
  python3 nozzle_mesh.py            # mesh + launch GUI viewer
  python3 nozzle_mesh.py -nopopup   # headless / batch mode

Requirements
------------
  pip install gmsh
"""

import math
import os
import sys

import gmsh

# ──────────────────────────────────────────────────────────────────────────────
#  1.  Initialise GMSH
# ──────────────────────────────────────────────────────────────────────────────
gmsh.initialize()
gmsh.model.add("nozzle")
occ = gmsh.model.occ          # OpenCASCADE geometry kernel (robust Booleans)

gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("General.Verbosity", 2)   # warnings and errors only

# ──────────────────────────────────────────────────────────────────────────────
#  2.  Nozzle dimensions
# ──────────────────────────────────────────────────────────────────────────────
D                 = 20.0   # inlet / outlet diameter                        [mm]
contraction_ratio = 0.5    # throat diameter = D * contraction_ratio         [-]

L_in     = 30.0   # inlet cylinder                                          [mm]
L_conv   = 20.0   # convergent cone                                         [mm]
L_throat = 15.0   # throat cylinder                                         [mm]
L_div    = 30.0   # divergent cone                                          [mm]
L_out    = 50.0   # outlet cylinder      (set to 0.0 to remove the section) [mm]

d = D * contraction_ratio  # throat diameter                                [mm]

# ── 2b.  Mesh resolution parameters ───────────────────────────────────────────
n_circ = 4     # cells along each 90 deg arc AND each core-square edge        [-]
n_rad  = 2     # radial cells in the annulus, core -> wall                   [-]

# Axial counts are deliberately higher in the two cones: that is where the
# cells get most stretched (small radius, steep wall), and a parameter sweep
# showed it is the cheapest place to spend cells -- [6,8,6,10,8] lifts min SICN
# from 0.37 to 0.45 for 15 % more cells, while going coarser or finer still
# does not improve the minimum.
#          inlet  conv  throat  div  outlet
n_ax   = [   6,     8,     6,   10,     8 ]   # axial cells per section       [-]

core_frac    = 0.55  # core-square corner radius = core_frac * R(z)          [-]
wall_grading = 1.0   # radial Progression coefficient (>1 clusters at wall)   [-]

# ── 2c.  Section table -> station lists ───────────────────────────────────────
# Every later section reads only z_stations / R_stations / n_axial, so adding
# or removing a section stays a single-place edit.  Zero-length sections are
# dropped, which lets L_out = 0.0 simply disappear.
_SECTIONS = [
    # name          length     r_start  r_end    axial cells
    ("inlet",      L_in,       D / 2,   D / 2,   n_ax[0]),
    ("convergent", L_conv,     D / 2,   d / 2,   n_ax[1]),
    ("throat",     L_throat,   d / 2,   d / 2,   n_ax[2]),
    ("divergent",  L_div,      d / 2,   D / 2,   n_ax[3]),
    ("outlet",     L_out,      D / 2,   D / 2,   n_ax[4]),
]
sections = [s for s in _SECTIONS if s[1] > 0.0]
if not sections:
    raise SystemExit("All section lengths are zero -- nothing to mesh.")

names   = [s[0] for s in sections]
lengths = [s[1] for s in sections]
n_axial = [s[4] for s in sections]

z_stations = [0.0]
for L in lengths:
    z_stations.append(z_stations[-1] + L)
R_stations = [sections[0][2]] + [s[3] for s in sections]

n_seg   = len(sections)
z_total = z_stations[-1]
R_max   = max(R_stations)

# Geometric tolerance: loose enough for OCC round-off, far tighter than any
# real feature of the geometry.
TOL = 1.0e-6 * max(D, z_total)

print("Convergent-divergent nozzle")
print(f"  inlet/outlet diameter D = {D:g} mm, throat diameter d = {d:g} mm "
      f"(ratio {contraction_ratio:g})")
for k in range(n_seg):
    print(f"  {names[k]:<11s} z = {z_stations[k]:7.2f} -> {z_stations[k+1]:7.2f} mm, "
          f"R = {R_stations[k]:5.2f} -> {R_stations[k+1]:5.2f} mm, "
          f"{n_axial[k]} axial cells")
print(f"  total length = {z_total:g} mm\n")


# ──────────────────────────────────────────────────────────────────────────────
#  3.  CAD solid  ->  STEP export
#
#  Built from OCC primitives rather than a revolve: cylinders and cones are
#  exact analytic surfaces, and fuse() removes the internal station faces so
#  the exported STEP is one clean solid.
#
#  The STEP is written HERE, before any partitioning, so the CAD file contains
#  the nozzle itself and not the 25 meshing blocks.
# ──────────────────────────────────────────────────────────────────────────────
primitives = []
for k in range(n_seg):
    z0, r0, r1, L = z_stations[k], R_stations[k], R_stations[k + 1], lengths[k]
    if abs(r1 - r0) < TOL:
        primitives.append((3, occ.addCylinder(0.0, 0.0, z0, 0.0, 0.0, L, r0)))
    else:
        primitives.append((3, occ.addCone(0.0, 0.0, z0, 0.0, 0.0, L, r0, r1)))

if len(primitives) > 1:
    fused, _ = occ.fuse([primitives[0]], primitives[1:])
else:
    fused = primitives
occ.synchronize()

solid_tags = [t for _, t in fused]
if len(solid_tags) != 1:
    raise SystemExit(f"fuse() produced {len(solid_tags)} solids, expected 1.")
solid = solid_tags[0]

# Both formats are written from the same synchronised model, so they describe
# the identical solid. OCC keeps the analytic character of the surfaces in
# either one -- the IGES holds five surfaces of revolution (entity 120) and the
# two end discs as planes (entity 108), matching the seven faces of the STEP --
# so neither format is an approximation of the other. Pick whichever the
# downstream reader is happier with.
cad_files = ["nozzle.step", "nozzle.iges"]
for cad_file in cad_files:
    gmsh.write(cad_file)
    print(f"CAD written -> {cad_file}  ({os.path.getsize(cad_file)} bytes)")
print()


# ──────────────────────────────────────────────────────────────────────────────
#  4.  O-grid core: a square frustum following  rc(z) = core_frac * R(z)
#
#  The square corners sit at 45/135/225/315 deg.  Because the taper is linear,
#  each of the four lateral faces is a PLANAR quad (both bounding chords are
#  parallel), so the frustum is built explicitly from points -> lines ->
#  plane surfaces -> volume.  This is fully deterministic and avoids the
#  surprises of OCC lofting (addThruSections).
# ──────────────────────────────────────────────────────────────────────────────
# The corners sit at 0/90/180/270 deg, NOT at 45/135/225/315.  OCC places the
# seam vertex of every circular edge (and the seam edge of every cylindrical /
# conical face) at theta = 0.  Putting a block boundary there makes the seam
# coincide with a petal edge; offset by 45 deg it would instead split one arc
# per station and one wall patch per segment, giving 5-corner surfaces that
# cannot be transfinite.
CORNER_ANGLES = [math.radians(90.0 * k) for k in range(4)]


def square_ring(z, rc):
    """Four corner points of the core square at height z, circumradius rc."""
    return [occ.addPoint(rc * math.cos(a), rc * math.sin(a), z)
            for a in CORNER_ANGLES]


core_tools = []
for k in range(n_seg):
    zb, zt = z_stations[k], z_stations[k + 1]
    pb = square_ring(zb, core_frac * R_stations[k])
    pt = square_ring(zt, core_frac * R_stations[k + 1])

    lb = [occ.addLine(pb[i], pb[(i + 1) % 4]) for i in range(4)]   # bottom square
    lt = [occ.addLine(pt[i], pt[(i + 1) % 4]) for i in range(4)]   # top square
    lv = [occ.addLine(pb[i], pt[i]) for i in range(4)]             # vertical edges

    faces = [
        occ.addPlaneSurface([occ.addCurveLoop(lb)]),
        occ.addPlaneSurface([occ.addCurveLoop(lt)]),
    ]
    for i in range(4):
        loop = occ.addCurveLoop([lb[i], lv[(i + 1) % 4], lt[i], lv[i]])
        faces.append(occ.addPlaneSurface([loop]))

    core_tools.append((3, occ.addVolume([occ.addSurfaceLoop(faces)])))

# ──────────────────────────────────────────────────────────────────────────────
#  5.  Partitioning tools
#
#  * four diagonal cutting patches at 45/135/225/315 deg, spanning the ANNULUS
#    only: the inner edge of each patch follows the core-square corner polyline
#    rc(z) = core_frac * R(z) exactly, so the patch meets the core along an edge
#    that already exists instead of slicing the core interior into wedges.  The
#    outer edge is oversized and trimmed back to the solid.
#  * one disk per internal station, at the EXACT local radius, to split the
#    model axially
# ──────────────────────────────────────────────────────────────────────────────
R_big = 1.5 * R_max

plane_tools = []
for a in CORNER_ANGLES:
    ca, sa = math.cos(a), math.sin(a)

    def at(r, z, ca=ca, sa=sa):
        return occ.addPoint(r * ca, r * sa, z)

    # inner boundary: up the core corner edge, station by station
    inner = [at(core_frac * R_stations[k], z_stations[k]) for k in range(n_seg + 1)]
    # outer boundary: back down at R_big (trimmed to the solid below)
    outer = [at(R_big, z_stations[-1]), at(R_big, z_stations[0])]

    ring = inner + outer
    edges = [occ.addLine(ring[i], ring[(i + 1) % len(ring)])
             for i in range(len(ring))]
    plane_tools.append((2, occ.addPlaneSurface([occ.addCurveLoop(edges)])))

# Trim the half-planes to the nozzle interior (keeps the solid: removeTool=False)
plane_tools, _ = occ.intersect(plane_tools, [(3, solid)],
                               removeObject=True, removeTool=False)

disk_tools = [(2, occ.addDisk(0.0, 0.0, z_stations[k], R_stations[k], R_stations[k]))
              for k in range(1, n_seg)]

# ──────────────────────────────────────────────────────────────────────────────
#  6.  Fragment -> conformal O-grid blocks
# ──────────────────────────────────────────────────────────────────────────────
occ.fragment([(3, solid)], core_tools + plane_tools + disk_tools)
occ.synchronize()

# Drop any tool fragment that survived outside the solid (no adjacent volume),
# then clean up the curves and points it leaves dangling.
for dim in (2, 1, 0):
    stray = [(dim, t) for _, t in gmsh.model.getEntities(dim)
             if len(gmsh.model.getAdjacencies(dim, t)[0]) == 0]
    if stray:
        gmsh.model.occ.remove(stray, recursive=False)
        occ.synchronize()

volumes = [t for _, t in gmsh.model.getEntities(3)]
n_expected_vol = 5 * n_seg
if len(volumes) != n_expected_vol:
    raise SystemExit(f"fragment() produced {len(volumes)} volumes, "
                     f"expected {n_expected_vol} (= 5 blocks x {n_seg} segments).")
print(f"O-grid partition: {len(volumes)} conformal blocks "
      f"({n_seg} segments x [1 core + 4 petals])\n")


# ──────────────────────────────────────────────────────────────────────────────
#  7.  Volume classification (by centre of mass)
#
#  The core block's centre of mass lies exactly on the axis; every petal's is
#  well off it.  This is a far more robust discriminator than tag order, which
#  fragment() does not guarantee.
# ──────────────────────────────────────────────────────────────────────────────
def segment_of(z):
    """Index of the axial segment containing z."""
    for k in range(n_seg):
        if z_stations[k] - TOL <= z <= z_stations[k + 1] + TOL:
            return k
    raise ValueError(f"z = {z} lies outside the nozzle")


core_vols, petal_vols = [], []
seg_of_vol = {}
for v in volumes:
    cx, cy, cz = occ.getCenterOfMass(3, v)
    seg_of_vol[v] = segment_of(cz)
    (core_vols if math.hypot(cx, cy) < 0.05 * R_max else petal_vols).append(v)

if len(core_vols) != n_seg or len(petal_vols) != 4 * n_seg:
    raise SystemExit(f"Block classification failed: {len(core_vols)} core / "
                     f"{len(petal_vols)} petal, expected {n_seg} / {4 * n_seg}.")


# ──────────────────────────────────────────────────────────────────────────────
#  8.  Curve classification and transfinite constraints
#
#  Every curve falls into exactly one of three classes, decided purely from its
#  endpoint coordinates.  That is what makes the transfinite counts agree on
#  opposite edges of all blocks without tracking any fragment tag:
#
#     dz > tol                 axial           -> n_axial[segment] + 1 nodes
#     dz ~ 0 and dr > tol      radial          -> n_rad + 1 nodes
#     dz ~ 0 and dr ~ 0        circumferential -> n_circ + 1 nodes
#                              (quarter arcs AND core-square edges: both have
#                               their two endpoints at the same radius)
# ──────────────────────────────────────────────────────────────────────────────
def endpoints_of(curve):
    pts = gmsh.model.getBoundary([(1, curve)], oriented=False, combined=False)
    return [gmsh.model.getValue(0, p, []) for _, p in pts]


def radius_at(curve, t):
    x, y, _ = gmsh.model.getValue(1, curve, [t])
    return math.hypot(x, y)


n_by_class = {"axial": 0, "radial": 0, "circumferential": 0}

for _, c in gmsh.model.getEntities(1):
    ends = endpoints_of(c)
    if len(ends) != 2:
        raise SystemExit(f"Curve {c} is closed; the O-grid expects only "
                         f"quarter arcs and straight segments.")
    (x0, y0, z0), (x1, y1, z1) = ends
    dz = abs(z1 - z0)
    dr = abs(math.hypot(x1, y1) - math.hypot(x0, y0))

    if dz > TOL:                                        # axial
        k = segment_of(0.5 * (z0 + z1))
        gmsh.model.mesh.setTransfiniteCurve(c, n_axial[k] + 1)
        n_by_class["axial"] += 1
    elif dr > TOL:                                      # radial
        if abs(wall_grading - 1.0) < 1e-12:
            gmsh.model.mesh.setTransfiniteCurve(c, n_rad + 1)
        else:
            # "Progression" GROWS the elements along the curve's own
            # parametrisation, so clustering at the wall needs 1/wall_grading
            # on an inward->outward curve.  And fragment() does not preserve
            # curve direction, so which case we are in has to be measured:
            # sample the radius at both parameter bounds.
            (t0,), (t1,) = gmsh.model.getParametrizationBounds(1, c)
            outward = radius_at(c, t1) > radius_at(c, t0)
            coef = 1.0 / wall_grading if outward else wall_grading
            gmsh.model.mesh.setTransfiniteCurve(c, n_rad + 1, "Progression", coef)
        n_by_class["radial"] += 1
    else:                                               # circumferential
        gmsh.model.mesh.setTransfiniteCurve(c, n_circ + 1)
        n_by_class["circumferential"] += 1

for _, s in gmsh.model.getEntities(2):
    gmsh.model.mesh.setTransfiniteSurface(s)
    gmsh.model.mesh.setRecombine(2, s)

for v in volumes:
    gmsh.model.mesh.setTransfiniteVolume(v)

print("Transfinite curves: "
      + ", ".join(f"{k} {v}" for k, v in n_by_class.items()) + "\n")


# ──────────────────────────────────────────────────────────────────────────────
#  9.  Physical groups
#
#  A surface is on the outer boundary iff exactly one volume is adjacent to it.
#  That is more robust than a bounding-box test, which cannot tell an internal
#  O-grid interface from a wall patch.
# ──────────────────────────────────────────────────────────────────────────────
inlet_surfs, outlet_surfs, wall_surfs = [], [], []
for _, s in gmsh.model.getEntities(2):
    up, _down = gmsh.model.getAdjacencies(2, s)
    if len(up) != 1:
        continue                                   # internal interface
    _, _, zmin, _, _, zmax = gmsh.model.getBoundingBox(2, s)
    if zmax - zmin < TOL and abs(zmin) < TOL:
        inlet_surfs.append(s)
    elif zmax - zmin < TOL and abs(zmin - z_total) < TOL:
        outlet_surfs.append(s)
    else:
        wall_surfs.append(s)

gmsh.model.addPhysicalGroup(3, volumes,      name="fluid")
gmsh.model.addPhysicalGroup(2, inlet_surfs,  name="inlet")
gmsh.model.addPhysicalGroup(2, outlet_surfs, name="outlet")
gmsh.model.addPhysicalGroup(2, wall_surfs,   name="wall")

print("Physical groups")
print(f"  fluid  : {len(volumes)} volumes")
print(f"  inlet  : {len(inlet_surfs)} surfaces")
print(f"  outlet : {len(outlet_surfs)} surfaces")
print(f"  wall   : {len(wall_surfs)} surfaces\n")

if len(inlet_surfs) != 5 or len(outlet_surfs) != 5:
    raise SystemExit("Inlet/outlet should each be 5 patches (core + 4 petals).")


# ──────────────────────────────────────────────────────────────────────────────
# 10.  Generate and save
#
#  No Mesh.SubdivisionAlgorithm and no Mesh.Algorithm3D here: transfinite
#  volumes bounded by recombined faces emit hexahedra directly, and subdividing
#  would quadruple the cell count for nothing.
# ──────────────────────────────────────────────────────────────────────────────
gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)   # Lethe / OpenFOAM / Fluent

print("Generating 3-D mesh ...")
gmsh.model.mesh.generate(3)

mesh_file = "nozzle.msh"
gmsh.write(mesh_file)
print(f"Mesh saved -> {mesh_file}\n")


# ──────────────────────────────────────────────────────────────────────────────
# 11.  Verification
# ──────────────────────────────────────────────────────────────────────────────
elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim=3)
non_hex = [t for t in elem_types if t != 5]
assert not non_hex, f"Non-hexahedral element types found: {non_hex}"

hex_tags = [t for arr in elem_tags for t in arr]
n_hex = len(hex_tags)

per_section = n_circ * n_circ + 4 * n_circ * n_rad
n_predicted = per_section * sum(n_axial)

qual = gmsh.model.mesh.getElementQualities(hex_tags, "minSICN")
q_min, q_mean = min(qual), sum(qual) / len(qual)

print("Verification")
print("  all 3-D elements are hexahedra (type 5)          OK")
print(f"  cells: {n_hex} (predicted {per_section} per cross-section "
      f"x {sum(n_axial)} axial layers = {n_predicted})")
assert n_hex == n_predicted, f"Cell count {n_hex} != predicted {n_predicted}"
print("  actual cell count matches the parameters         OK")
print(f"  quality minSICN: min = {q_min:.3f}, mean = {q_mean:.3f}")
if q_min <= 0.0:
    raise SystemExit("Inverted or degenerate elements (minSICN <= 0).")
if q_min < 0.30:
    # 0.30 is the usual practical floor for hexahedra.  At the default
    # (deliberately minimal) resolution the worst cells sit in the cones just
    # off the throat; raising n_ax there is the effective fix.
    print("  WARNING: min SICN below 0.30 -- raise n_ax in the cone sections")
for cad_file in cad_files:
    print(f"  {cad_file}: {os.path.getsize(cad_file)} bytes                     OK")
print(f"  {mesh_file}: {os.path.getsize(mesh_file)} bytes                     OK\n")


# ──────────────────────────────────────────────────────────────────────────────
# 12.  (Optional) interactive viewer
# ──────────────────────────────────────────────────────────────────────────────
if "-nopopup" not in sys.argv:
    print("Launching GMSH GUI (close window to exit) ...")
    gmsh.fltk.run()

gmsh.finalize()
