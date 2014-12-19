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
 * @defgroup boundary Boundary and manifold description for triangulations
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
 * object has an invalid manifold id. On the other hand, the boundary
 * id of faces and the manifold id of objects can be set either at
 * creation time or later by looping over all cells and querying their
 * faces.
 *
 * It is then possible to associate objects describing the geometry to
 * certain boundary_id values used in a triangulation and to certain
 * manifold_id values.
 *
 * Befor version 8.2, the library allowed only boundary faces to
 * follow a curved geometric description. Since version 8.2 this has
 * been introduced also for interior faces and cells, and the
 * boundary_id has been separated from the manifold_id.
 *
 * Although the old behavior is still supported, one should use the
 * boundary indicator only for the physical meaning associated, for
 * example, to boundary conditions, and revert to manifold_ids to
 * describe the geometry of the triangulation.
 *
 * The behavior of the Triangulation class w.r.t. geometry
 * descriptions is the following: Triangulation::set_boundary() and
 * Triangulation::set_manifold() do the exact same thing: they attach
 * a manifold descriptor to the specified id. The first function
 * expects a Boundary descriptor (which is a specialization of a
 * Manifold description) and is provided mainly for backward
 * compatibility, while the second class expects a Manifold
 * descriptor. Notice that the Triangulation class only uses the
 * Manifold interface, and you could describe both the interior and
 * the boundary of the domain using the same object. The additional
 * information contained in the Boundary interface is related to the
 * computation of the exact normals. 
 *
 * Whenever a new vertex is needed in an object, the Triangulation
 * queries the manifold_id of the object which needs refinement. If
 * the manifold_id is set to numbers::invalid_manifold_id, then the
 * Triangulation queries the boundary_id (if the face is on the
 * boundary) or the material_id (if the Triangulation is of
 * codimension one and the object is a cell). If the previous queries
 * resulted in a number different from numbers::invalid_manifold_id,
 * then the Triangulation looks wether a previous call to
 * Triangulation::set_manifold() (or set_boundary()) was performed
 * with the given id, and if yes, it uses the stored object to obtain
 * new vertices, otherwise it uses a FlatManifold or StraightBoundary
 * object.
 *
 * @note This behavior is backward compatible to that of deal.II versions
 * prior to 8.2. If one ignores the manifold_id of an object (i.e., if it has
 * never been set), by default it is and remains set to
 * numbers::invalid_manifold_id. In that case, the first query above will
 * trigger a query to the old style boundary_id. This behavior will be
 * maintained for a while, but might eventually be changed. The suggested
 * strategy is to use manifold_ids to describe the geometry, and boundary_ids
 * to describe boundary conditions.
 * 
 *
 * @see @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
 * @see @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, Luca Heltai, 1998-2014
 */
