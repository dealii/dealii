// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.templates.h>


DEAL_II_NAMESPACE_OPEN

#include "lac/vector_memory.inst"
template class VectorMemory<
  LinearAlgebra::distributed::Vector<float, MemorySpace::Default>>;
template class VectorMemory<
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>;
template class GrowingVectorMemory<
  LinearAlgebra::distributed::Vector<float, MemorySpace::Default>>;
template class GrowingVectorMemory<
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>;

namespace internal
{
  namespace GrowingVectorMemoryImplementation
  {
    void
    release_all_unused_memory()
    {
#include "lac/vector_memory_release.inst"
      dealii::GrowingVectorMemory<dealii::LinearAlgebra::distributed::Vector<
        float,
        MemorySpace::Default>>::release_unused_memory();
      dealii::GrowingVectorMemory<dealii::LinearAlgebra::distributed::Vector<
        double,
        MemorySpace::Default>>::release_unused_memory();
    }
  } // namespace GrowingVectorMemoryImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
