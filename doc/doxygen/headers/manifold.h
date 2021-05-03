// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


/**
 * @defgroup manifold Manifold description for triangulations
 *
 * <h3>Overview</h3>
 *
 * The classes in this module are concerned with the description of the
 * manifold in which the domain that a Triangulation describes lives. This
 * manifold description is necessary in several contexts:
 *
 * <ul>
 *
 *   <li> Mesh refinement: Whenever a cell is refined, it is necessary
 *   to introduce new vertices in the Triangulation. In the
 *   simplest case, one assumes that the objects that make up the
 *   Triangulation are straight line segments, a bi-linear surface or
 *   a tri-linear volume. The next vertex is then simply put into the
 *   middle of the old ones (where "middle" means a suitable average of the
 *   locations of the pre-existing vertices). This is the default behavior of
 *   the Triangulation class, and is described by the FlatManifold class.
 *
 *   On the other hand, if one deals with curved geometries, or geometries which
 *   require a denser refinement in some direction, this is not the appropriate
 *   thing to do. The classes derived from the Manifold base class therefore
 *   describe the geometry of a domain. One can then attach an object of a class
 *   derived from this base class to the Triangulation object using the
 *   Triangulation::set_manifold() function associating it with a manifold_id
 *   (see types::manifold_id), use this manifold_id on the cells, faces or edges
 *   of the triangulation that should be described by this manifold using the
 *   TriaAccessor::set_manifold_id() function, and then the Triangulation will
 *   ask the manifold object where a new vertex to be located on a cell, face or
 *   edge so attributed should be located upon mesh refinement. Several classes
 *   already exist to support the most common geometries, e.g.,
 *   CylindricalManifold, or PolarManifold, which represent respectively the
 *   geometry obtained when describing your space in cylindrical coordinates or
 *   in polar coordinates. By default, all curved geometries generated using
 *   functions in the GridGenerator namespace attach the correct Manifold object
 *   to the curved parts of the domain.
 *
 *   <li> Integration: When using higher order finite element methods, it is
 *   often necessary to compute cell terms (like cell contributions to the
 *   matrix and right hand side of the linear system) using curved
 *   approximations of the boundary, rather than the straight line
 *   approximation. The actual implementation of such curved elements happens
 *   in the Mapping class (see the @ref mapping module), which however obtains
 *   its information about the boundary of the domain from the classes
 *   described here. The same is, of course, true when integrating boundary
 *   terms (e.g., inhomogeneous Neumann boundary conditions).
 *
 *   <li> Domains with nonzero codimension: In cases where a Triangulation is
 *   embedded into a higher dimensional space, i.e., whenever the second
 *   template argument of the Triangulation class is explicitly specified and
 *   larger than the first (for an example, see step-34), the manifold
 *   description objects serve as a tool to describe the geometry not only of
 *   the boundary of the domain but of the domain itself, in case the domain
 *   is a manifold that is in fact curved. In these cases, one can use the
 *   Triangulation::set_manifold() function to indicate what manifold
 *   description to use when refining the curve, or when computing integrals
 *   using high order mappings.
 *
 * </ul>
 * Many other examples, as well as much theoretical underpinning for the
 * implementation in deal.II, is provided in the
 * @ref geometry_paper "geometry paper".
 *
 * In deal.II, a Manifold is seen as a collection of points, together
 * with a notion of distance between points (on the manifold). New
 * points are typically obtained by providing a local coordinate
 * system on the manifold, identifying existing points in the local
 * coordinate system (pulling them back using the local map to obtain
 * their local coordinates), find the new point in the local
 * coordinate system by weighted sums of the existing points, and
 * transforming back the point in the real space (pushing it forward
 * using the local map). The main class that implements this mechanism
 * is the ChartManifold class, and this is the class that users will
 * likely overload for complex geometries.
 *
 * While this process is non trivial in most cases of interest, for most of
 * the trivial geometries, like cylinders, spheres or shells, deal.II provides
 * reasonable implementations. More complicated examples can be described
 * using the techniques shown in step-53 and step-54.
 *
 * In the grand scheme of things, the classes of this module interact
 * with a variety of other parts of the library:
 * @dot
 digraph G
{
  graph[rankdir="TB",bgcolor="transparent"];

  node [fontname="FreeSans",fontsize=15,
        shape=box,height=0.2,width=0.4,
        color="black", fillcolor="white", style="filled"];
  edge [color="black", weight=10];

  tria       [label="Triangulation",    URL="\ref grid"];
  fe         [label="Finite elements",    URL="\ref feall"];
  mapping    [label="Mapping",          URL="\ref mapping"];
  quadrature [label="Quadrature",       URL="\ref Quadrature"];
  dh         [label="DoFHandler",       URL="\ref dofs"];
  fevalues   [label="FEValues",         URL="\ref feaccess"];
  systems    [label="Linear systems",   URL="\ref LAC"];
  solvers    [label="Linear solvers",   URL="\ref Solvers"];
  output     [label="Graphical output", URL="\ref output"];
  manifold   [label="Manifold",         URL="\ref manifold", fillcolor="deepskyblue"];

  tria -> dh              [color="black",style="solid"];
  fe -> dh                [color="black",style="solid"];
  fe -> fevalues          [color="black",style="solid"];
  mapping -> fevalues     [color="black",style="solid"];
  quadrature -> fevalues  [color="black",style="solid"];
  dh -> systems           [color="black",style="solid"];
  fevalues -> systems     [color="black",style="solid"];
  systems -> solvers      [color="black",style="solid"];
  solvers -> output       [color="black",style="solid"];
  manifold -> tria        [color="black",style="solid"];
  manifold -> mapping     [color="black",style="solid"];

  {
    rank=same
    mapping -> quadrature [dir="none", color="transparent"];
    quadrature -> fe      [dir="none", color="transparent"];
    fe -> tria            [dir="none", color="transparent"];
  }

  node [fontname="FreeSans",fontsize=12,
        shape=record,height=0.2,width=0.4,
        color="gray55", fontcolor="gray55", fillcolor="white", style="filled"];
  edge [color="gray55", weight=1];

  opencascade [label="OpenCASCADE"];
  opencascade -> manifold [dir="none"];


  node [fontname="FreeSans",fontsize=12,
        shape=ellipse,height=0.2,width=0.4,
        color="gray55", fontcolor="gray55", fillcolor="white", style="filled"];
  edge [color="gray55", weight=1];

  gmsh        [label="gmsh", URL="\ref Gmsh"];
  gmsh -> tria       [dir="none"];
}
 * @enddot
 *
 *
 * <h3>An example</h3>
 *
 * A simple example why dealing with curved geometries is already provided by
 * step-1, though it is not elaborated there. By default, the functions in
 * GridGenerator will attach manifolds to meshes when needed. In each code
 * snippet below we call Triangulation::reset_all_manifolds() to remove these
 * manifolds and handle all Manifold attachment in the example itself to make
 * the impact of the choice of Manifold clear.
 *
 * Consider this
 * small variation of the <code>second_grid()</code> function shown there,
 * where we simply refine <i>every</i> cell several times:
 * @code
 *  const Point<2> center (1,0);
 *  const double inner_radius = 0.5,
 *               outer_radius = 1.0;
 *  Triangulation<2> triangulation;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              10);
 *  // as noted above: disable all non-Cartesian manifolds
 *  // for demonstration purposes:
 *  triangulation.reset_all_manifolds();
 *
 *  triangulation.refine_global (3);
 * @endcode
 * This code leads to a mesh that looks like this:
 *
 * @image html hypershell-nothing.png ""
 *
 * Our intention was to get a mesh that resembles a ring. However, since we did
 * not describe this to the triangulation, what happens is that we start with
 * the 10 coarse cells in circumferential direction we told
 * GridGenerator::hyper_shell() to create, and each of these is then 3 times
 * globally refined. Each time refinement requires a new vertex, it is placed
 * in the middle of the existing ones, regardless of what we may have intended
 * (but omitted to describe in code).
 *
 * This is easily remedied. Consider this code:
 * @code
 *  const Point<2> center (1,0);
 *  const SphericalManifold<2> manifold(center);
 *  const double inner_radius = 0.5,
 *               outer_radius = 1.0;
 *  Triangulation<2> triangulation;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              10);
 *  // again disable all manifolds for demonstration purposes
 *  triangulation.reset_all_manifolds();
 *  triangulation.set_all_manifold_ids_on_boundary(0);
 *  triangulation.set_manifold (0, manifold);
 *
 *  triangulation.refine_global (3);
 * @endcode
 * This code is better, producing the following mesh:
 *
 * @image html hypershell-boundary-only.png ""
 *
 * The mesh looks better in that it faithfully reproduces the circular inner
 * and outer boundaries of the domain. However, it is still possible to
 * identify 20 kinks in the tangential lines. They result from the fact that
 * every time a cell is refined, new vertices on interior lines are just
 * placed into the middle of the existing line (the boundary lines are handled
 * differently because we have attached a manifold object). In the first
 * refinement with 10 cells, we got improved points because both outer
 * boundaries have provided a curved description according to the description
 * on blending different manifolds below. In other words, the new points after
 * the first refinement end up in places that may be in the geometric middle
 * of a straight line, but not on a circle around the center.
 *
 * This can be remedied by assigning a manifold description not only to
 * the lines along the boundary, but also to the radial lines and cells (which,
 * in turn, will inherit it to the new lines that are created upon mesh
 * refinement). This is exactly what GridGenerator::hyper_shell() does by default.
 * For demonstration purposes, we disable the default Manifold behavior and then
 * duplicate it manually:
 * @code
 *  const Point<2> center (1,0);
 *  const SphericalManifold<2> manifold(center);
 *  const double inner_radius = 0.5,
 *               outer_radius = 1.0;
 *  Triangulation<2> triangulation;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              10);
 *  // again disable all manifolds for demonstration purposes
 *  triangulation.reset_all_manifolds();
 *  // reenable the manifold:
 *  triangulation.set_all_manifold_ids(0);
 *  triangulation.set_manifold (0, manifold);
 *  triangulation.refine_global (3);
 * @endcode
 * This leads to the following mesh:
 *
 * @image html hypershell-all.png ""
 *
 * So why does this matter? After all, the last two meshes describe the
 * exact same domain and we know that upon mesh refinement we obtain the
 * correct solution regardless of the choice of cells, as long as the
 * diameter of the largest cell goes to zero.
 *
 * There are two answers to this question. First, the numerical effort
 * of solving a partial differential equation to a certain accuracy typically
 * depends on the <i>quality</i> of cells since the constant $C$ in error
 * estimates of the form $\|u-u_h\|_{H^1} \le Ch^p \|u\|_{H^{p+1}}$ depends
 * on factors such as the maximal ratio of radii of the smallest circumscribed
 * to largest inscribed circle over all cells (for triangles; or a suitable
 * generalization for other types of cells). Thus, it is worthwhile creating
 * meshes with cells that are as well-formed as possible. This is arguably
 * not so much of an issue for the meshes shown above, but is sometimes an
 * issue. Consider, for example, the following code and mesh:
 * @code
 *  const Point<2> center (1,0);
 *  const SphericalManifold<2> manifold(center);
 *  const double inner_radius = 0.5,
 *               outer_radius = 1.0;
 *  Triangulation<2> triangulation;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              3);    // three circumferential cells
 *  triangulation.reset_all_manifolds();
 *  triangulation.set_all_manifold_ids_on_boundary(0);
 *  triangulation.set_manifold (0, manifold);
 *
 *  triangulation.refine_global (3);
 * @endcode
 *
 * @image html hypershell-boundary-only-3.png ""
 *
 * Here, we create only three circumferential cells in the beginning, and
 * refining them leads to the mesh shown. Clearly, we have cells with bad
 * aspect ratios, despite the first refinement that puts the new point into
 * the middle.
 *
 * If we drive this further and start with a coarse mesh of a much thinner rim
 * between the radii 0.8 and 1.0 and only three cells (which
 * may be inappropriate here, since we know that it is not sufficient, but may
 * also be impossible to avoid for complex geometries generated in mesh
 * generators), we observe the following:
 *
 * @code
 *  const Point<2> center (1,0);
 *  const SphericalManifold<2> manifold(center);
 *  const double inner_radius = 0.8,
 *               outer_radius = 1.0;
 *  Triangulation<2> triangulation;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              3);    // three circumferential cells
 *  triangulation.reset_all_manifolds();
 *  triangulation.set_all_manifold_ids_on_boundary(0);
 *  triangulation.set_manifold (0, manifold);
 *
 *  triangulation.refine_global (3);
 * @endcode
 *
 * @image html hypershell-boundary-thin-3.png ""
 *
 * This mesh neither has the correct geometry after refinement, nor do
 * all cells have positive area as is necessary for the finite element
 * method to work. However, even when starting with such an inopportune
 * mesh, we can make things work by attaching a suitable geometry description
 * not only to the boundary but also to interior cells and edges, using
 * the same code as above:
 * @code
 *  const Point<2> center (1,0);
 *  const double inner_radius = 0.8,
 *               outer_radius = 1.0;
 *  Triangulation<2> triangulation;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              3);    // three circumferential cells
 *
 *  triangulation.refine_global (3);
 * @endcode
 *
 * @image html hypershell-all-3.png ""
 *
 * In this last example we finally let GridGenerator do its job and we keep
 * the default manifold configuration, which is a SphericalManifold on every
 * cell and face.
 *
 * Here, even starting with an initial, inappropriately chosen mesh retains
 * our ability to adequately refine the mesh into one that will serve us
 * well. This example may be manufactured here, but it is relevant, for example
 * in the context of what GridGenerator::hyper_shell() produces in 3d
 * (see the documentation of this function). It is also germane to the
 * cases discussed in the @ref GlossDistorted "glossary entry on distorted cells".
 *
 * @see @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
 *
 * <h3>Computing the weights for combining different manifold descriptions</h3>
 *
 * In a realistic application, it happens regularly that different manifold
 * descriptions need to be combined. The simplest case is when a curved
 * description is only available for the boundary but not for the interior of
 * the computational domain. The manifold description for a ball also falls
 * into this category, as it needs to combine a spherical manifold at the
 * circular part with a straight-sided description in the center of the domain
 * where the spherical manifold is not valid.
 *
 * In general, the process of blending different manifold descriptions in
 * deal.II is achieved by the so-called transfinite interpolation. Its formula
 * in 2D is, for example, described on <a
 * href="https://en.wikipedia.org/wiki/Transfinite_interpolation">
 * Wikipedia</a>. Given a point $(u,v)$ on a chart, the image of this point
 * in real space is given by
 * @f{align*}{
 * \mathbf S(u,v) &= (1-v)\mathbf c_0(u)+v \mathbf c_1(u) + (1-u)\mathbf c_2(v) + u \mathbf c_3(v) \\
 * &\quad - \left[(1-u)(1-v) \mathbf x_0 + u(1-v) \mathbf x_1 + (1-u)v \mathbf x_2 + uv \mathbf x_3 \right]
 * @f}
 * where $\bf x_0, \bf x_1, \bf x_2, \bf x_3$ denote the four vertices
 * bounding the image space and $\bf c_0, \bf c_1, \bf c_2, \bf c_3$ are the
 * four curves describing the lines of the cell.
 *
 * If we want to find the center of the cell according to the manifold (that
 * is also used when the grid is refined), the chart is the unit cell
 * $(0,1)^2$ and we want to evaluate this formula in the point $(u,v) = (0.5,
 * 0.5)$. In that case, $\mathbf c_0(0.5)$ is the position of the midpoint of
 * the lower face (indexed by 2 in deal.II's ordering) that is derived from
 * its own manifold, $\mathbf c_1(0.5)$ is the position of the midpoint of the
 * upper face (indexed by 3 in deal.II), $\mathbf c_2(0.5)$ is the midpoint of
 * the face on the left (indexed by 0), and $\mathbf c_3(0.5)$ is the midpoint
 * of the right face. In this formula, the weights equate to
 * $\frac{\displaystyle 1}{\displaystyle 2}$ for the four midpoints in the
 * faces and to $-\frac{\displaystyle 1}{\displaystyle 4}$ for the four
 * vertices. These weights look weird at first sight because the vertices
 * enter with negative weight but the mechanism does what we want: In case of
 * a cell with curved description on two opposite faces but straight lines on
 * the other two faces, the negative weights of $-\frac{\displaystyle
 * 1}{\displaystyle 4}$ in the vertices balance with the center of the two
 * straight lines in radial direction that get weight $\frac{\displaystyle
 * 1}{\displaystyle 2}$. Thus, the average is taken over the two center points
 * in curved direction, exactly placing the new point in the middle.
 *
 * In three spatial dimensions, the weights are $+\frac{\displaystyle
 * 1}{\displaystyle 2}$ for the face midpoints, $-\frac{\displaystyle
 * 1}{\displaystyle 4}$ for the line mid points, and $\frac{\displaystyle
 * 1}{\displaystyle 8}$ for the vertices, again balancing the different
 * entities. In case all the surrounding of a cell is straight, the formula
 * reduces to the obvious weight $\frac{\displaystyle 1}{\displaystyle 8}$ on
 * each of the eight vertices.
 *
 * In the MappingQGeneric class, a generalization of this concept to the
 * support points of a polynomial representation of curved cells, the nodes of
 * the Gauss-Lobatto quadrature, is implemented by evaluating the boundary
 * curves in the respective Gauss-Lobatto points $(u_i,v_i)$ and combining
 * them with the above formula. The weights have been verified to yield
 * optimal convergence rates $\mathcal O(h^{k+1})$ also for very high
 * polynomial degrees, say $k=10$.
 *
 * In the literature, other boundary descriptions are also used. Before
 * version 9.0 deal.II used something called Laplace smoothing where the
 * weights that are applied to the nodes on the circumference to get the
 * position of the interior nodes are determined by solving a Laplace equation
 * on the unit element. However, this led to boundary layers close to the
 * curved description, i.e., singularities in the higher derivatives of the
 * mapping from unit to real cell.
 *
 * If the transition from a curved boundary description to a straight
 * description in the interior is done wrong, it is typically impossible to
 * achieve high order convergence rates. For example, the Laplace smoothing
 * inside a single cell leads to a singularity in the fourth derivative of the
 * mapping from the reference to the real cell, limiting the convergence rate
 * to 3 in the cells at the boundary (and 3.5 if global L2 errors were
 * measured in 2D). Other more crude strategies, like completely ignoring the
 * presence of two different manifolds and simply computing the additional
 * points of a high-order mapping in a straight coordinate system, could lead
 * to even worse convergence rates. The current implementation in deal.II, on
 * the other hand, has been extensively verified in this respect and should
 * behave optimally.
 *
 * A bad strategy for blending a curved boundary representation with flat
 * interior representations obviously also reflects mesh quality. For example,
 * the above case with only 3 circumferential cells leads to the following
 * mesh with Laplace manifold smoothing rather than the interpolation from the
 * boundary as is implemented in deal.II:
 *
 * @image html hypershell-boundary-only-3-old.png ""
 *
 * To use a more practical example, consider the refinement of a ball with a
 * SphericalManifold attached to the spherical surface. The Laplace-type smoothing
 * gives the following rather poor mesh:
 *
 * @image html hyperball-mesh-smoothing-laplace.png ""
 *
 * If we, instead, use the weights derived from transfinite interpolation, the
 * situation is considerably improved:
 *
 * @image html hyperball-mesh-smoothing-interpolate.png ""
 *
 * Of course, one could get even better meshes by applying the
 * TransfiniteInterpolationManifold to the whole domain except the boundary
 * where SphericalManifold is attached, as shown by the figures in that class,
 * but in principle, the mesh smoothing implemented in deal.II is as good as
 * it can get from a boundary description alone.
 *
 * @ingroup grid
 * @author Luca Heltai, 2013, Martin Kronbichler, 2017
 */
