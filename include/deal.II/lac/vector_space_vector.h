// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2023 by the deal.II authors
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
