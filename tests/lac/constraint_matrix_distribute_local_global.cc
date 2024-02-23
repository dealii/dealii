// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <stdio.h>

#include <iostream>

#include "../tests.h"

// Tests distribute_local_to_global for
// (1) the case that a local matrix
//     with all diagonals equal to zero is
//     distributed while a dof is
//     constrained, and
// (2) the case that all entries of the local matrix
//     are zero while a dof is constrained.


int
main()
{
  initlog();

  // set up constraint
  AffineConstraints<double> constraints;
  constraints.constrain_dof_to_zero(0);
  constraints.close();

  // global matrix and vector
  FullMatrix<double> global_matrix(2);
  Vector<double>     global_vector(2);

  // first test: add matrix with all diagonals zero
  FullMatrix<double> local_matrix(2);
  local_matrix(1, 0) = local_matrix(0, 1) = 1.;
  Vector<double> local_vector(2);
  constraints.distribute_local_to_global(
    local_matrix, local_vector, {0, 1}, global_matrix, global_vector, true);
  // output result
  for (unsigned int m = 0; m < global_matrix.m(); ++m)
    {
      for (unsigned int n = 0; n < global_matrix.n(); ++n)
        deallog << global_matrix(m, n) << ' ';
      deallog << std::endl;
    }
  deallog << std::endl;

  // second test: add matrix with all entries zero
  global_matrix = 0.;
  local_matrix  = 0.;
  constraints.distribute_local_to_global(
    local_matrix, local_vector, {0, 1}, global_matrix, global_vector, true);
  // output result
  for (unsigned int m = 0; m < global_matrix.m(); ++m)
    {
      for (unsigned int n = 0; n < global_matrix.n(); ++n)
        deallog << global_matrix(m, n) << ' ';
      deallog << std::endl;
    }
  deallog << std::endl;
}
