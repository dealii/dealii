// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.templates.h>


DEAL_II_NAMESPACE_OPEN

template class VectorMemory<LinearAlgebra::CUDAWrappers::Vector<float>>;
template class VectorMemory<LinearAlgebra::CUDAWrappers::Vector<double>>;
template class GrowingVectorMemory<LinearAlgebra::CUDAWrappers::Vector<float>>;
template class GrowingVectorMemory<LinearAlgebra::CUDAWrappers::Vector<double>>;

namespace internal
{
  namespace GrowingVectorMemoryImplementation
  {
    void
    release_all_unused_cuda_memory()
    {
      dealii::GrowingVectorMemory<dealii::LinearAlgebra::CUDAWrappers::Vector<
        float>>::release_unused_memory();
      dealii::GrowingVectorMemory<dealii::LinearAlgebra::CUDAWrappers::Vector<
        double>>::release_unused_memory();
    }
  } // namespace GrowingVectorMemoryImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
