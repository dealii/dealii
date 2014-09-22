// ---------------------------------------------------------------------
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
 * @defgroup hp hp finite element support
 *
 * Classes and functions that have to do with hp finite elements. The step-27
 * tutorial program gives an overview of how to use the classes in this
 * namespace. A slightly more exotic application is given in step-46.
 *
 * The hp namespace implements the algorithms and data structures used for
 * the hp framework in deal.II. An overview over the details of how these
 * algorithms work and what data structures are used is given in the
 * @ref hp_paper "hp paper".
 */

/**
 * @defgroup hpcollection hp Collections
 *
 * In the implementation of the hp finite element method, each cell might have
 * a different finite element associated with it. To handle this, the
 * hp::DoFHandler must have a whole set of finite element classes associated
 * with it. This concept is represented by the hp::FECollection class: Objects
 * of this type act as containers that hold a whole set of finite element
 * objects. Instead of storing pointers to finite element objects on each
 * cell, we then only store an index for each cell that identifies the finite
 * element object within the collection that should be used by this cell. The
 * DoFHandler object associated with the given cell can then assign degrees of
 * freedom to each cell in accordance with the finite element used for it.
 *
 * A similar situation arises when integrating terms on a cell: one may want
 * to use different quadrature formulas for different finite elements. For
 * example, on cells where we use a Q1 element, a QGauss(2) object (i.e. a
 * quadrature formula with two points in each space direction) may be
 * sufficient, but on another cell where a Q3 element is used, this would lead
 * to underintegration and we should use a QGauss(4) formula instead. Just as
 * above, there exists a class hp::QCollection that acts as a collection of
 * quadrature formulas
 *
 * Finally, one may want to use different orders for the boundary
 * approximation for cells with different orders for the finite element. The
 * hp::MappingCollection class allows to do this.
 *
 * All of these three classes, the hp::FECollection, hp::QCollection,
 * and hp::MappingCollection classes, implement an interface very
 * similar to that of <code>std::vector</code>. They have functions
 * <code>push_back()</code> to add a finite element, quadrature
 * formula, or mapping to the collection. They have an
 * <code>operator[] (unsigned int)</code> function that allows to
 * retrieve a reference to a given element of the collection. And they
 * have a <code>size()</code> function that returns the number of
 * elements in the collection. Some of the classes, in particular that
 * holding finite element objects, also implement other functions
 * specific to their purpose.
 *
 * The similarity goes beyond the interface: When adding an element to the
 * collection, all of the classes create a copy of the argument. This allows
 * to pass a temporary object to the function adding the element. For example,
 * the following works:
 * @verbatim
 *   FECollection<dim> fe_collection;
 *   for (unsigned int degree=1; degree<5; ++degree)
 *     fe_collection.push_back (FE_Q<dim>(degree));
 * @endverbatim
 * 
 * This way, one can add elements of polynomial degree 1 through 4 to the
 * collection. It is not necessary to retain the added object: the collection
 * makes a coyp of it, it does not only store a pointer to the given finite
 * element object. This same observation also holds for the other collection
 * classes.
 *
 * It is customary that within an hp finite element program, one keeps
 * collections of finite elements and quadrature formulas with the same number
 * of elements, each element of the one collection matching the element in the
 * other. This is not necessary, but it often makes coding a lot simpler. If a
 * collection of mappings is used, the same holds for hp::MappingCollection
 * objects as well.
 *
 * @ingroup hp
 */


/**
 * A namespace for the implementation of hp finite element specific algorithms
 * and data structures.
 *
 * @ingroup hp
 */
namespace hp
{
}
