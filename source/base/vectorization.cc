// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/vectorization.h>

DEAL_II_NAMESPACE_OPEN

// VectorizedArray must be a POD (plain old data) type to make sure it
// can use maximum level of compiler optimization.
// A type is POD if it has standard layout (similar to a C struct)
// and it is trivial (can be statically default initialized)
// Here, the trait std::is_pod cannot be used because it is deprecated
// in C++20.
//
// Check these statements to ensure we catch problems if we accidentally
// make these classes non-POD.
static_assert(
  std::is_standard_layout_v<VectorizedArray<double>> &&
    std::is_trivially_default_constructible_v<VectorizedArray<double>> &&
    std::is_trivially_copyable_v<VectorizedArray<double>>,
  "VectorizedArray<double> must be a POD type");
static_assert(
  std::is_standard_layout_v<VectorizedArray<float>> &&
    std::is_trivially_default_constructible_v<VectorizedArray<float>> &&
    std::is_trivially_copyable_v<VectorizedArray<float>>,
  "VectorizedArray<float> must be a POD type");

DEAL_II_NAMESPACE_CLOSE
