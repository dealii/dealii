// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

#include <deal.II/matrix_free/cuda_matrix_free.templates.h>

#ifdef DEAL_II_WITH_CUDA

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  // Do not instantiate for dim = 1
  template class MatrixFree<2, float>;
  template class MatrixFree<2, double>;
  template class MatrixFree<3, float>;
  template class MatrixFree<3, double>;
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
