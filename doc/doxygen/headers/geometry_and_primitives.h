// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
 * In deal.II, Triangulation objects are built from line segments, triangles or
 * quadrilaterals, or tetrahedra, pyramids, wedges, or hexahedra (depending on
 * the space dimension). The ReferenceCell class encodes all properties of these
 * basic objects in (as the name implies) reference coordinates, such as number
 * of vertices per cell, lines per face, and the coordinates of each vertex.
 * This abstraction enables writing applications mostly independently of the
 * actual space dimension as well as the ReferenceCell types of a Triangulation:
 * i.e., if you consistently use ReferenceCell's data members, then the same
 * program should work with both quadrilateral and triangular cells. For
 * example, loops over all cell vertices would simply run from zero to
 * `cell->reference_cell().n_vertices()` instead of hard-coding values which
 * only work with a particular type of reference cell in a particular dimension.
 * In this way, the program will be correct for all reference cell types, and
 * one can run a program in a different space dimension simply by recompilation
 * instead of having to change a significant portion of the code. These
 * dimension-independent programming techniques are extensively discussed in the
 * first few tutorial programs and are used throughout deal.II.
 */
