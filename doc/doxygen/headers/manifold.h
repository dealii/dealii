// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


/**
 * @defgroup manifold Manifold description for triangulations
 *
 * <h3>Overview</h3>
 *
 * The classes in this module are concerned with the description of the
 * manifold in which the domain that a Triangulation describes lives. This
 * manifold description is necessary in two contexts:
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
 *   On the other hand, if one deals with curved geometries, or geometries
 *   which require a denser refinement in some direction, this is not the
 *   appropriate thing to do. The classes derived from the Manifold base class
 *   therefore describe the geometry of a domain. One can then attach an
 *   object of a class derived from this base class to the Triangulation
 *   object using the Triangulation::set_manifold() function associating it
 *   with a manifold id (see types::manifold_id), use this manifold id on the
 *   cells, faces or edges of the triangulation that should be described by
 *   this manifold using the TriaAccessor::set_manifold_id() function, and
 *   then the Triangulation will ask the manifold object where a new vertex to
 *   be located on a cell, face or edge so attributed should be located upon
 *   mesh refinement. Several classes already exist to support the most common
 *   geometries, e.g., CylinderManifold, or PolarManifold, which represent
 *   respectively the geometry obtained when describing your space in
 *   cylindrical coordintes or in polar coordinates.
 *
 *   <li> Integration: When using higher order finite element methods, it is
 *   often necessary to compute cell terms (like cell contributions to the
 *   matrix and right hand side of the linear system) using curved
 *   approximations of the boundary, rather than the straight line
 *   approximation. The actual implementation of such curved elements happens
 *   in the Mapping class (see the @ref mapping module), which however obtains
 *   its information about the boundary of the domain from the classes
 *   described here. The same is, of course, true when integrating boundary
 *   terms (e.g., inhomogenous Neumann boundary conditions).
 *
 * </ul>
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
 * using the techniques shown in step-53.
 *
 * The boundary of a Triangulation is a special case of Manifold, for
 * which additional information can be useful in user codes, such as
 * normal vectors to surfaces or to curves. If your coarse mesh is reasonably
 * shaped, you might be interested in only attaching a manifold
 * description to boundary portion of your domain. This can be done
 * using the Triangulation::set_boundary() function, which take as arguments a Boundary
 * object (derived from Manifold). Notice that Triangulation uses only
 * the Manifold interface, not the Boundary interface. Other tools,
 * however, might need to compute exact normals at quadrature points,
 * and therefore a wrapper to query Boundary objects is provided. 
 *
 *
 * <h3>An example</h3>
 *
 * A simple example why dealing with curved geometries is already provided
 * by step-1, though it is not elaborated there. For example, consider this
 * small variation of the <code>second_grid()</code> function shown there,
 * where we simply refine <i>every</i> cell several times and do not deal
 * with boundaries at all:
 * @code
 *  Triangulation<2> triangulation;
 *
 *  const Point<2> center (1,0);
 *  const double inner_radius = 0.5,
 *               outer_radius = 1.0;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              10);
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
 * This is easily remedied. step-1 already shows how to do this. Consider this
 * code:
 * @code
 *  Triangulation<2> triangulation;
 *
 *  const Point<2> center (1,0);
 *  const double inner_radius = 0.5,
 *               outer_radius = 1.0;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              10);
 *  const HyperShellBoundary<2> boundary_description(center);
 *  triangulation.set_boundary (0, boundary_description);
 *
 *  triangulation.refine_global (3); 
 * @endcode
 * This code is better, producing the following mesh:
 *
 * @image html hypershell-boundary-only.png ""
 *
 * The mesh looks better in that it faithfully reproduces the circular inner
 * and outer boundaries of the domain. However, it is still possible to
 * identify the original 10 cells by the kinks in the tangential lines. They
 * result from the fact that every time a cell is refined, new vertices on
 * interior lines are just placed into the middle of the existing line (the
 * boundary lines are handled differently because we have attached boundary
 * objects). In other words, they end up in places that may be in the geometric
 * middle of a straight line, but not on a circle around the center.
 *
 * This can be remedied by assigning a manifold description not only to
 * the lines along the boundary, but also to the radial lines and cells (which,
 * in turn, will inherit it to the new lines that are created upon mesh
 * refinement). This code achieves this:
 * @code
 *  Triangulation<2> triangulation;
 *
 *  const Point<2> center (1,0);
 *  const double inner_radius = 0.5,
 *               outer_radius = 1.0;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              10);
 *  const SphericalManifold<2> boundary_description(center);
 *  triangulation.set_manifold (0, boundary_description);
 *
 *  Triangulation<2>::active_cell_iterator
 *    cell = triangulation.begin_active(),
 *    endc = triangulation.end();
 *  for (; cell!=endc; ++cell)
 *    cell->set_all_manifold_ids (0);  
 *  
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
 *  Triangulation<2> triangulation;
 *
 *  const Point<2> center (1,0);
 *  const double inner_radius = 0.5,
 *               outer_radius = 1.0;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              4);    // four circumferential cells
 *  const HyperShellBoundary<2> boundary_description(center);
 *  triangulation.set_boundary (0, boundary_description);
 *
 *  triangulation.refine_global (3); 
 * @endcode
 *
 * @image html hypershell-boundary-only-4.png ""
 *
 * Here, we create only four circumferential cells in the beginning,
 * and refining them leads to the mesh shown. Clearly, here we have
 * cells with bad aspect ratios.
 *
 * If we drive this further and start with a coarse mesh of only
 * three cells (which may be inappropriate here, since we know that it
 * is not sufficient, but may also be impossible to avoid for complex
 * geometries generated in mesh generators), then we obtain the following
 * mesh:
 *
 * @image html hypershell-boundary-only-3.png ""
 *
 * This mesh neither has the correct geometry after refinement, nor do
 * all cells have positive area as is necessary for the finite element
 * method to work. However, even when starting with such in inopportune
 * mesh, we can make things work by attaching a suitable geometry description
 * not only to the boundary but also to interior cells and edges, using
 * the same code as above:
 * @code
 *  Triangulation<2> triangulation;
 *
 *  const Point<2> center (1,0);
 *  const double inner_radius = 0.5,
 *               outer_radius = 1.0;
 *  GridGenerator::hyper_shell (triangulation,
 *                              center, inner_radius, outer_radius,
 *                              3);    // three circumferential cells
 *  const SphericalManifold<2> boundary_description(center);
 *  triangulation.set_manifold (0, boundary_description);
 *
 *  Triangulation<2>::active_cell_iterator
 *    cell = triangulation.begin_active(),
 *    endc = triangulation.end();
 *  for (; cell!=endc; ++cell)
 *    cell->set_all_manifold_ids (0);  
 *  
 *  triangulation.refine_global (3);
 * @endcode
 *
 * @image html hypershell-all-3.png ""
 *
 * Here, even starting with an initial, inappropriately chosen mesh retains
 * our ability to adequately refine the mesh into one that will serve us
 * well. This example may be manifactured here, but it is relevant, for example
 * in the context of what GridGenerator::hyper_shell() produces in 3d
 * (see the documentation of this function). It is also germane to the
 * cases discussed in the @ref GlossDistorted "glossary entry on distorted cells".
 * 
 *
 * @see @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
 *
 * @ingroup grid
 * @author Luca Heltai, 2013
 */
