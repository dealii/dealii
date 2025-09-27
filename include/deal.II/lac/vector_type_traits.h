// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
 *   struct is_serial_vector<VectorType> : std::true_type
 *   {};
 * @endcode
 * for a serial vector type, respectively,
 * @code
 *   template <>
 *   struct is_serial_vector<VectorType> : std::false_type
 *   {};
 * @endcode
 * for a vector type with support of distributed storage,
 * must be done in a header file of a vector declaration.
 */
template <typename T>
struct is_serial_vector;


DEAL_II_NAMESPACE_CLOSE

#endif
