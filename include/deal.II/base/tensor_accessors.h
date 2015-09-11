// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

#ifndef dealii__tensor_accessors_h
#define dealii__tensor_accessors_h

#include <deal.II/base/config.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/table_indices.h>


DEAL_II_NAMESPACE_OPEN

/**
 * This namespace is a collection of algorithms working on generic
 * tensorial objects (of arbitrary rank).
 *
 * The rationale to implement such functionality in a generic fashion in a
 * separate namespace is
 *  - to easy code reusability and therefore avoid code duplication.
 *  - to have a well-defined interface that allows to exchange the low
 *    level implementation.
 *
 *
 * A tensorial object has the notion of a rank and allows a rank-times
 * recursive application of the index operator, e.g., if <code>t</code> is
 * a tensorial object of rank 4, the following access is valid:
 * @code
 *   t[1][2][1][4]
 * @endcode
 *
 * deal.II has its own implementation for tensorial objects such as
 * dealii::Tensor<rank, dim, Number> and
 * dealii::SymmetricTensor<rank, dim, Number>
 *
 * The methods and algorithms implemented in this namespace, however, are
 * fully generic. More precisely, it can operate on nested c-style arrays,
 * or on class types <code>T</code> with a minimal interface that provides
 * a local typedef <code>value_type</code> and an index operator
 * <code>operator[](unsigned int)</code> that returns a (const or
 * non-const) reference of <code>value_type</code>:
 * @code
 *   template<...>
 *   class T
 *   {
 *     typedef ... value_type;
 *     value_type & operator[](unsigned int);
 *     const value_type & operator[](unsigned int) const;
 *   };
 * @endcode
 *
 * This namespace provides primitves for access, reordering and contraction
 * of such objects.
 *
 * @ingroup geomprimitives
 */
namespace TensorAccessors
{


} /* namespace TensorAccessors */

DEAL_II_NAMESPACE_CLOSE

#endif /* dealii__tensor_accessors_h */
