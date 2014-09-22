// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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
 * @defgroup geomprimitives Geometric and other primitives
 *
 * This group contains a number of classes that act as geometric
 * primitives or primitives for other mathematical objects. For
 * example, the Tensor @<rank,dim@> class provides tensors of rank
 * <code>rank</code> in <code>dim</code> space dimensions. Likewise,
 * the SymmetricTensor offers symmetric tensors.
 *
 * Geometrically, the Point class is the foundation of all geometric
 * descriptions in the deal.II library. It denotes a geometric point
 * in <code>dim</code> dimensional space. One can view a point as a
 * vector the with <code>dim</code> coordinates that connects the
 * origin with that particular point; as such, the Point class is
 * derived from tensors of rank 1 (i.e. vectors), but in contrast to
 * arbitrary tensors points have the special connotation of points in
 * space, and therefore have some additional properties.
 *
 * In deal.II, meshes are built from line segments, quadrilaterals, or
 * hexahedra (depending on the space dimension). The GeometryInfo
 * class is used to describe properties of these basic objects in unit
 * space (i.e. for the unit line, unit square, and unit cube). It
 * offers static data members denoting the number of vertices per
 * cell, lines per face, or where which vertex is located. This
 * abstraction allows to write applications mostly independently of
 * the actual space dimension: loops over all vertices would simply
 * run from zero to GeometryInfo<dim>::vertices_per_cell instead of
 * from 0 to 4 (in 2d) or 0 to 8 (in 3d). In this way, the program
 * will be correct in 2d as well as 3d, and one can run a program in a
 * different space dimension simply by recompilation instead of having
 * to change a significant portion of the code. These
 * dimension-independent programming techniques are extensively
 * discussed in the first few tutorial programs and are used
 * throughout deal.II.
 */

