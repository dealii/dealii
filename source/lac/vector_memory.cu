// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.templates.h>


DEAL_II_NAMESPACE_OPEN

template class VectorMemory<LinearAlgebra::CUDAWrappers::Vector<float>>;
template class VectorMemory<LinearAlgebra::CUDAWrappers::Vector<double>>;
template class GrowingVectorMemory<LinearAlgebra::CUDAWrappers::Vector<float>>;
template class GrowingVectorMemory<LinearAlgebra::CUDAWrappers::Vector<double>>;
template class VectorMemory<
  LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA>>;
template class VectorMemory<
  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>>;
template class GrowingVectorMemory<
  LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA>>;
template class GrowingVectorMemory<
  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>>;

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
      dealii::GrowingVectorMemory<
        dealii::LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA>>::
        release_unused_memory();
      dealii::GrowingVectorMemory<
        dealii::LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>>::
        release_unused_memory();
    }
  } // namespace GrowingVectorMemoryImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
