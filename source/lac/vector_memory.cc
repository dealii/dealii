// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2019 by the deal.II authors
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

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.templates.h>


DEAL_II_NAMESPACE_OPEN

#include "vector_memory.inst"

namespace internal
{
  namespace GrowingVectorMemoryImplementation
  {
#ifdef DEAL_II_WITH_CUDA
    void
    release_all_unused_cuda_memory();
#endif

    void
    release_all_unused_memory()
    {
#include "vector_memory_release.inst"
#ifdef DEAL_II_WITH_CUDA
      release_all_unused_cuda_memory();
#endif
    }
  } // namespace GrowingVectorMemoryImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
