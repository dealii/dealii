// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2023 by the deal.II authors
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

#ifndef dealii_cuda_fe_evaluation_h
#define dealii_cuda_fe_evaluation_h

#include <deal.II/matrix_free/portable_fe_evaluation.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace for the CUDA wrappers
 */
// GCC 9 and before do not recognize the [[deprecated]] attribute
#if defined(__GNUC__) && (__GNUC__ < 10)
namespace CUDAWrappers
#else
namespace DEAL_II_DEPRECATED_EARLY CUDAWrappers
#endif
{
  using namespace Portable;
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
