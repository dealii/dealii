// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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
 * @defgroup boundary Boundary description for triangulations
 *
 * The classes in this module are concerned with the description of the
 * boundary of a domain in which a Triangulation lives. This boundary
 * description is necessary in three contexts:
 * <ul>
 *
 *   <li> Mesh refinement: Whenever a cell at the boundary is refined, it is
 *   necessary to introduce at least one new vertex on the boundary. In the
 *   simplest case, one assumes that the boundary consists of straight line
 *   segments (in 2d) or a bilinear surface (in 3d) between the vertices of
 *   the original, coarsest mesh, and the next vertex is simply put into the
 *   middle of the old ones. This is the default behavior of the Triangulation
 *   class, and is described by the StraightBoundary class.
 *
 *   On the other hand, if one deals with curved boundaries, this is not the
 *   appropriate thing to do. The classes derived from the Boundary base class
 *   therefore describe the geometry of a domain. One can then attach an
 *   object of a class derived from this base class to the Triangulation
 *   object using the Triangulation::set_boundary() function, and the
 *   Triangulation will ask the boundary object where a new vertex should be
 *   located upon mesh refinement. Several classes already exist to support
 *   the most common geometries, e.g., CylinderBoundary, HyperBallBoundary, or
 *   HyperShellBoundary.
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
 *   <li> In cases where a Triangulation is embedded into a higher dimensional
 *   space, i.e., whenever the second template argument of the Triangulation
 *   class is explicitly specified and larger than the first (for an example,
 *   see step-34), the boundary description objects also serve as a tool to
 *   describe the geometry not only of the boundary of the domain but of the
 *   domain itself, in case the domain is a manifold that is in fact curved.
 * </ul>
 *
 * In the context of triangulations, each face of a cell that is located at
 * the boundary of the domain stores a number called <tt>boundary_id</tt> that
 * uniquely identifies which part of the boundary this face is on. If nothing
 * is specified at creation time, each boundary face has a zero boundary
 * id. On the other hand, the boundary id of faces can be set either at
 * creation time or later by looping over all cells and querying their faces.
 *
 * It is then possible to associate objects describing the boundary to certain
 * boundary_id values used in a triangulation. Note that this is not
 * necessary: in some cases one may want to use the default straight boundary
 * approximation, and use non-zero boundary indicators for completely
 * different purposes, for example to indicate that a part of the boundary has
 * a different kind of boundary condition in the partial differential
 * equation.
 *
 * @see @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, 1998-2006
 */
