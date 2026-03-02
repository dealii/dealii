// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.templates.h>
#include <deal.II/lac/lapack_full_matrix.h>

DEAL_II_NAMESPACE_OPEN

#include "lac/la_parallel_block_vector.inst"

// do a few functions that currently don't fit the scheme because they have
// two template arguments that need to be different (the case of same
// arguments is covered by the default copy constructor and copy operator that
// is declared separately)

namespace LinearAlgebra
{
  namespace distributed
  {
#ifndef DOXYGEN
#  define TEMPL_COPY_CONSTRUCTOR(S1, S2)                    \
    template BlockVector<S1, ::dealii::MemorySpace::Host> & \
    BlockVector<S1, ::dealii::MemorySpace::Host>::operator= \
      <S2>(const BlockVector<S2, ::dealii::MemorySpace::Host> &)

    TEMPL_COPY_CONSTRUCTOR(double, float);
    TEMPL_COPY_CONSTRUCTOR(float, double);

#  undef TEMPL_COPY_CONSTRUCTOR

    template void
    BlockVector<float, ::dealii::MemorySpace::Default>::reinit<float>(
      const BlockVector<float, ::dealii::MemorySpace::Default> &,
      const bool);
    template void
    BlockVector<double, ::dealii::MemorySpace::Default>::reinit<double>(
      const BlockVector<double, ::dealii::MemorySpace::Default> &,
      const bool);
#endif
  } // namespace distributed
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE
