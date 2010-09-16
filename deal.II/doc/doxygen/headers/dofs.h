//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @defgroup dofs Degrees of Freedom
 *
 * This module groups classes and namespaces that have to do with
 * handling degrees of freedom. The central class of this group is the
 * DoFHandler class: it is built on top of a triangulation and a
 * finite element class and allocated degrees of freedom on each cell
 * of the triangulation as required for the finite element space
 * described by the finite element object. There are other variants of
 * the DoFHandler class such as hp::DoFHandler,
 * parallel::distributed::DoFHandler and MGDoFHandler that do similar
 * things for more special cases.
 *
 * DoFHandlers extend Triangulation objects (and the other classes in the @ref
 * grid module) in that they, too, offer iterators that run over all cells,
 * faces, or other geometric objects that make up a triangulation. These
 * iterators are derived from the triangulation iterators and therefore offer
 * the same functionality, but they also offer additional functions. For
 * example, they allow to query the indices of the degrees of freedom
 * associated with the present cell.
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
