// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#include <deal.II/base/vectorization.h>

DEAL_II_NAMESPACE_OPEN

// VectorizedArray must be a POD (plain old data) type to make sure it
// can use maximum level of compiler optimization.
// A type is POD if it has standard layout (similar to a C struct)
// and it is trivial (can be statically default initialized)
// Here, the trait std::is_pod cannot be used because it is deprecated
// in C++20.
static_assert(std::is_standard_layout<VectorizedArray<double>>::value &&
                std::is_trivial<VectorizedArray<double>>::value,
              "VectorizedArray<double> must be a POD type");
static_assert(std::is_standard_layout<VectorizedArray<float>>::value &&
                std::is_trivial<VectorizedArray<float>>::value,
              "VectorizedArray<float> must be a POD type");

// For the specializations of VectorizedArray, we need to instantiate the
// static constexpr variable for some compilers. On the other hand, MSCV wants
// us explicitly to not do so, otherwise we get: "error C2908: explicit
// specialization 'const unsigned int
// dealii::VectorizedArray<double, 2>::n_array_elements' has already been
// instantiated".
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && !defined(DEAL_II_MSVC)
#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
const unsigned int VectorizedArray<double, 8>::n_array_elements;
const unsigned int VectorizedArray<float, 16>::n_array_elements;
#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
const unsigned int VectorizedArray<double, 4>::n_array_elements;
const unsigned int VectorizedArray<float, 8>::n_array_elements;
#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
const unsigned int VectorizedArray<double, 2>::n_array_elements;
const unsigned int VectorizedArray<float, 4>::n_array_elements;
#  endif
#endif

DEAL_II_NAMESPACE_CLOSE
