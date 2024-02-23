// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_space_vector_h
#define dealii_vector_space_vector_h

#include <deal.II/base/config.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN

// ---------------------------- Free functions --------------------------

namespace LinearAlgebra
{
  /**
   * Shift all entries of the vector by a constant factor so that the mean
   * value of the vector becomes zero.
   */
  template <typename VectorType>
  void
  set_zero_mean_value(VectorType &vector)
  {
    vector.add(-vector.mean_value());
  }
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
