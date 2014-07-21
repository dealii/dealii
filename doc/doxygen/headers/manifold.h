// ---------------------------------------------------------------------
// $Id$
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
 * The classes in this module are concerned with the description of the
 * manifold in which the domain that a Triangulation describes lives. This
 * manifold description is necessary in two contexts:
 *
 * <ul>
 *
 *   <li> Mesh refinement: Whenever a cell is refined, it is necessary
 *   to introduce some new vertices in the Triangulation. In the
 *   simplest case, one assumes that the objects that make up the
 *   Triangulation are straight line segments, a bi-linear surface or
 *   a tri-linear volume. The next vertex is then simply put into the
 *   middle of the old ones. This is the default behavior of the
 *   Triangulation class, and is described by the FlatManifold class.
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
 * In deal.II, a Manifold is seen as a collection of points, together with a
 * notion of distance between points (on the manifold). New points are
 * typically obtained by providing a local coordinate system on the manifold,
 * identifying existing points in the local coordinate system (pulling them
 * back using the local map to obtain their local coordinates), find the new
 * point in the local coordinate system by weighted sums of the existing
 * points, and transforming back the point in the real space (pushing it
 * forward using the local map).  (However, this pull back and push forward
 * only happens in classes that implement manifolds, and your own class can
 * also choose other ways of figuring out where a new point should be
 * located.)
 *
 * While this process is non trivial in most cases of interest, for
 * most of the trivial geometries, like cylinders, spheres or shells,
 * we provide reasonable implementations. 
 *
 * @see @ref GlossManifoldIndicator "Glossary entry on manifold
 * indicators"
 *
 * @ingroup grid
 * @author Luca Heltai, 2013
 */
