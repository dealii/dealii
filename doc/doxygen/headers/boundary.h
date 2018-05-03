// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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
 * @defgroup boundary Boundary and manifold description for triangulations
 *
 * @warning This module describes usage of the Boundary classes, which are
 * deprecated: see the
 * @ref manifold
 * module for detailed information on how to describe curved boundaries with
 * the Manifold classes, which are the replacements for the Boundary classes.
 *
 * The classes in this module are concerned with the description of the
 * geometry of a domain in which a Triangulation lives. This geometry
 * description is necessary in three contexts:
 * <ul>
 *
 *   <li> Mesh refinement: Whenever a cell is refined, it is necessary
 *   to introduce at least one new vertex. In the simplest case, one
 *   assumes that the cells and their faces consists of straight line
 *   segments, bilinear surface or trilinear volumes between the
 *   vertices of the original, coarsest mesh, and the next vertex is
 *   simply put into the middle of the old ones. This is the default
 *   behavior of the Triangulation class, and is described by the
 *   StraightBoundary and FlatManifold classes.
 *
 *   On the other hand, if one deals with curved geometries and
 *   boundaries, this is not the appropriate thing to do. The classes
 *   derived from the Manifold and Boundary base classes describe the
 *   geometry of a domain. One can then attach an object of a class
 *   derived from this base classes to the Triangulation object using
 *   the Triangulation::set_boundary() or
 *   Triangulation::set_manifold() functions, and the Triangulation
 *   will ask the manifold object where a new vertex should be located
 *   upon mesh refinement. Several classes already exist to support
 *   the boundary of the most common geometries, e.g.,
 *   CylinderBoundary, HyperBallBoundary, or HyperShellBoundary.
 *
 *   <li> Integration: When using higher order finite element methods,
 *   it is often necessary to compute cell terms (like cell
 *   contributions to the matrix and right hand side of the linear
 *   system) using curved approximations of the geometry, rather than
 *   the straight line approximation. The actual implementation of
 *   such curved elements happens in the Mapping class (see the @ref
 *   mapping module), which however obtains its information about the
 *   manifold description from the classes described here. The same
 *   is, of course, true when integrating boundary terms (e.g.,
 *   inhomogenous Neumann boundary conditions).
 *
 *   <li> In cases where a Triangulation is embedded into a higher
 *   dimensional space, i.e., whenever the second template argument of
 *   the Triangulation class is explicitly specified and larger than
 *   the first (for an example, see step-34), the manifold description
 *   objects serve as a tool to describe the geometry not only of the
 *   boundary of the domain but of the domain itself, in case the
 *   domain is a manifold that is in fact curved. In these cases, one
 *   can use the Triangulation::set_manifold() function to indicate
 *   what manifold description to use when refining the curve, or when
 *   computing integrals using high order mappings.
 *
 * </ul>
 *
 * In the context of triangulations, each object stores a number
 * called <tt>manifold_id</tt>, and each face of a cell that is
 * located at the boundary of the domain stores a number called
 * <tt>boundary_id</tt> that uniquely identifies which part of the
 * boundary this face is on. If nothing is specified at creation time,
 * each boundary face has a zero boundary id and each triangulation
 * object has a flat_manifold_id. On the other hand, the boundary
 * id of faces and the manifold_id of objects can be set either at
 * creation time or later by looping over all cells and querying their
 * faces.
 *
 * It is then possible to associate objects describing the geometry to certain
 * manifold_id values.
 *
 * Before version 8.2, the library allowed only boundary faces to follow a
 * curved geometric description. Since version 8.2 this has been introduced also
 * for interior faces and cells, and the boundary_id has been separated from the
 * manifold_id. The former is used to identify the type of boundary conditions
 * to apply, while the latter is used to identify the geometry and describe how
 * new vertices should be created upon refinement of the mesh, or where high
 * order Mapping objects should place their support points on the exact
 * geometry.
 *
 * Since version 9.0 of the library, the boundary_id associated to the boundary
 * faces is ignored by Manifold objects, and Manifold descriptors can only be
 * attached to manifold ids.
 *
 * The behavior of the Triangulation class with respect to geometry descriptions
 * is the following: Triangulation::set_manifold() attaches a manifold
 * descriptor to the specified manifold_id. The function expects a Manifold
 * descriptor, and you could describe both the interior and the boundary of the
 * domain using the same object.
 *
 * Whenever a new vertex is needed in an object, the Triangulation queries the
 * manifold_id of the object which needs refinement. If the query resulted in a
 * number different from numbers::flat_manifold_id, then the Triangulation looks
 * whether a previous call to Triangulation::set_manifold() was performed with
 * the given id. If it was, then the triangulation uses the stored object to
 * obtain new vertices, otherwise it uses a FlatManifold object.
 *
 * @note This behavior is **not** backward compatible to that of deal.II
 * versions prior to 9.0. If one ignores the manifold_id of an object (i.e., if
 * it has never been set), by default it is and remains set to
 * numbers::flat_manifold_id. In previous versions of the library, the
 * boundary_id of the boundary faces would be queried and used in place of the
 * manifold_id. If you have old programs that only set boundary ids, you should
 * modify them to use manifold ids instead (or you could use
 * GridTools::copy_boundary_to_manifold_ids or
 * GridTools::map_boundary_to_manifold_ids)
 *
 * @see @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
 * @see @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, Luca Heltai, 1998-2018
 */
