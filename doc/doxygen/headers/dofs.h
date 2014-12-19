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
 * @defgroup dofs Degrees of Freedom
 *
 * This module groups classes and namespaces that have to do with
 * handling degrees of freedom. The central class of this group is the
 * DoFHandler class: it is built on top of a triangulation and a
 * finite element class and allocated degrees of freedom on each cell
 * of the triangulation as required for the finite element space
 * described by the finite element object. There are other variants of
 * the DoFHandler class such as hp::DoFHandler that do similar
 * things for more special cases.
 *
 * DoFHandler objects are used together with objects of type FiniteElement
 * (or hp::FECollection in the case of hp::DoFHandler) to enumerate all the
 * degrees of freedom that exist in a triangulation for this particular
 * finite element. As such, the combination of mesh, finite element, and
 * DoF handler object can be thought of as providing a <i>basis</i> of
 * the finite element space: the mesh provides the locations at which basis
 * functions are defined; the finite element describes what kinds of basis
 * functions exist; and the DoF handler object provides an enumeration of
 * the basis, i.e., it is provides a concrete structure of the space so that
 * we can describe functions in this finite dimensional space by vectors
 * of coefficients.
 *
 * DoFHandlers extend Triangulation objects (and the other classes in the @ref
 * grid module) in that they, too, offer iterators that run over all cells,
 * faces, or other geometric objects that make up a triangulation. These
 * iterators are derived from the triangulation iterators and therefore offer
 * the same functionality, but they also offer additional functions. For
 * example, they allow to query the indices of the degrees of freedom
 * associated with the present cell. Note that DoFHandler classes are <i>not
 * derived</i> from Triangulation, though they use Triangulation objects;
 * the reason is that there can be more than one DoFHandler object that works
 * on the same Triangulation object.
 *
 * In addition to the DoF handler classes, this module holds a number of
 * auxiliary classes not commonly used in application programs, as well as
 * three classes that are not directly associated with the data structures of
 * the DoFHandler class. The first of these is the ConstraintMatrix class that
 * stores and treats the constraints associated with hanging nodes. Secondly,
 * the DoFRenumbering namespace offers functions that can reorder degrees of
 * freedom; among its functions are ones that sort degrees of freedom in
 * downstream direction, for example, and ones that sort degrees of freedom in
 * such a way that the bandwidth of associated matrices is minimized. Finally,
 * the DoFTools class offers a variety of algorithms around handling degrees
 * of freedom.
 */
