// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2023 by the deal.II authors
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
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_base.h>


DEAL_II_NAMESPACE_OPEN


template <typename VectorType>
void
MGSmootherBase<VectorType>::apply(const unsigned int level,
                                  VectorType        &u,
                                  const VectorType  &rhs) const
{
  u = typename VectorType::value_type(0.);
  smooth(level, u, rhs);
}



template <typename VectorType>
void
MGTransferBase<VectorType>::prolongate_and_add(const unsigned int to_level,
                                               VectorType        &dst,
                                               const VectorType  &src) const
{
  VectorType temp;
  temp.reinit(dst, true);

  this->prolongate(to_level, temp, src);

  dst += temp;
}


// Explicit instantiations

#include "mg_base.inst"

DEAL_II_NAMESPACE_CLOSE
