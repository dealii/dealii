// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/lac/vector_memory.templates.h>

DEAL_II_NAMESPACE_OPEN

#include "vector_memory.inst"

#ifdef DEAL_II_WITH_CUDA
template class VectorMemory<LinearAlgebra::CUDAWrappers::Vector<float>>;
template class VectorMemory<LinearAlgebra::CUDAWrappers::Vector<double>>;
template class GrowingVectorMemory<LinearAlgebra::CUDAWrappers::Vector<float>>;
template class GrowingVectorMemory<LinearAlgebra::CUDAWrappers::Vector<double>>;
#endif

namespace internal
{
  namespace GrowingVectorMemoryImplementation
  {
    void
    release_all_unused_memory()
    {
#include "vector_memory_release.inst"

#ifdef DEAL_II_WITH_CUDA
      dealii::GrowingVectorMemory<dealii::LinearAlgebra::CUDAWrappers::Vector<
        float>>::release_unused_memory();
      dealii::GrowingVectorMemory<dealii::LinearAlgebra::CUDAWrappers::Vector<
        double>>::release_unused_memory();
#endif
    }
  } // namespace GrowingVectorMemoryImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
