// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_vector_type_traits_h
#define dealii_vector_type_traits_h

#include <deal.II/base/config.h>

#include <type_traits>


DEAL_II_NAMESPACE_OPEN


/**
 * Type trait for a serial vector, i.e. a vector class for which storage is not
 * supported to be distributed over processes.
 *
 * The specialization
 * @code
 *   template <>
 *   struct is_serial_vector< VectorType > : std::true_type {};
 * @endcode
 * for a serial vector type, respectively,
 * @code
 *   template <>
 *   struct is_serial_vector< VectorType > : std::false_type {};
 * @endcode
 * for a vector type with support of distributed storage,
 * must be done in a header file of a vector declaration.
 *
 * @author Uwe Koecher, 2017
 */
template <typename T>
struct is_serial_vector;


DEAL_II_NAMESPACE_CLOSE

#endif
