
// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
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


  for (unsigned int row = 0; row < ghost_penalty_1d.n_rows(); ++row)
    {
      for (unsigned int col = 0; col < ghost_penalty_1d.n_cols(); ++col)
        deallog << std::setw(10) << std::fixed << std::setprecision(4)
                << ghost_penalty_1d(row, col) << " ";

      deallog << std::endl;
    }
  deallog << std::endl;
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