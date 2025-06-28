
// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// check creation of 1d ghost penalty matrix


#include <deal.II/numerics/tensor_product_matrix_creator.h>

#include "../tests.h"

void
test(const unsigned int fe_degree)
{
  FullMatrix<double> ghost_penalty_1d =
    TensorProductMatrixCreator::create_1d_ghost_penalty_matrix(
      FE_Q<1>(fe_degree), 1.);


  std::stringstream sstring;

  sstring << "Ghost penalty matrix for FE_Q<1>(" << fe_degree << "):\n";
  ghost_penalty_1d.print_formatted(sstring, 4, true, 10, "0.");

  deallog << sstring.str() << std::endl;

  // Check that the matrix is symmetric
  for (unsigned int i = 0; i < ghost_penalty_1d.m(); ++i)
    for (unsigned int j = 0; j < ghost_penalty_1d.n(); ++j)
      AssertThrow(ghost_penalty_1d(i, j) == ghost_penalty_1d(j, i),
                  ExcMessage("Matrix is not symmetric"));
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  test(1);
  test(3);
}
