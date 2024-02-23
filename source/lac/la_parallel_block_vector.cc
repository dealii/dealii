// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.templates.h>
#include <deal.II/lac/lapack_full_matrix.h>

DEAL_II_NAMESPACE_OPEN


// disable instantiation for MSVC for now because of a compiler bug,
// see https://github.com/dealii/dealii/issues/2875
#ifndef DEAL_II_MSVC

#  include "la_parallel_block_vector.inst"

// do a few functions that currently don't fit the scheme because they have
// two template arguments that need to be different (the case of same
// arguments is covered by the default copy constructor and copy operator that
// is declared separately)

namespace LinearAlgebra
{
  namespace distributed
  {
#  define TEMPL_COPY_CONSTRUCTOR(S1, S2)                 \
    template BlockVector<S1> &BlockVector<S1>::operator= \
      <S2>(const BlockVector<S2> &)

    TEMPL_COPY_CONSTRUCTOR(double, float);
    TEMPL_COPY_CONSTRUCTOR(float, double);

#  undef TEMPL_COPY_CONSTRUCTOR
  } // namespace distributed
} // namespace LinearAlgebra

#endif // ! DEAL_II_MSVC

DEAL_II_NAMESPACE_CLOSE
