
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

// chceck creation of 1D mass and laplace matrices.

#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/numerics/tensor_product_matrix_creator.h>

#include "../tests.h"


void
print_matrix(const FullMatrix<double> &matrix, const std::string &label)
{
  std::stringstream sstring;
  sstring << label << "\n";
  matrix.print_formatted(sstring, 4, true, 10, "0.");
  deallog << sstring.str() << std::endl;
}

template <int fe_degree>
void
test()
{
  FE_Q<1>   fe_q(fe_degree);
  FE_DGQ<1> fe_dgq(fe_degree);
  auto      mass_matrix =
    TensorProductMatrixCreator::create_1d_cell_mass_matrix(fe_q,
                                                           1.,
                                                           {true, true});
  print_matrix(mass_matrix,
               "Mass matrix for FE_Q<1>(" + std::to_string(fe_degree) + "):");



  FEEvaluation<1, fe_degree, fe_degree + 1, 1, double> fe_eval(
    fe_q, QGauss<1>(fe_degree + 1), update_values);


  auto mass_matrix_renumbered =
    TensorProductMatrixCreator::create_1d_cell_mass_matrix(
      fe_q, 1., {true, true}, fe_eval.get_internal_dof_numbering());
  print_matrix(mass_matrix_renumbered,
               "Mass matrix for FE_Q<1>(" + std::to_string(fe_degree) +
                 ") renumbered with FEEvaluation:");

  auto mass_matrix_dgq =
    TensorProductMatrixCreator::create_1d_cell_mass_matrix(fe_dgq,
                                                           1.,
                                                           {true, true});

  print_matrix(mass_matrix_dgq,
               "Mass matrix for FE_DGQ<1>(" + std::to_string(fe_degree) + "):");

  for (unsigned int i = 0; i < mass_matrix_dgq.m(); ++i)
    for (unsigned int j = 0; j < mass_matrix_dgq.n(); ++j)
      if (std::abs(mass_matrix_dgq(i, j) - mass_matrix_renumbered(i, j)) >
          1e-12)
        deallog << "Different at " << i << " " << j << std::endl;



  auto mass_matrix_patch =
    TensorProductMatrixCreator::create_1D_discretization_matrix(
      mass_matrix_renumbered, 2, 1, {true, true});

  print_matrix(mass_matrix_patch, "Mass matrix patch:");


  auto laplace_matrix =
    TensorProductMatrixCreator::create_1d_cell_laplace_matrix(fe_q,
                                                              1.,
                                                              {true, true});
  print_matrix(laplace_matrix,
               "Laplace matrix for FE_Q<1>(" + std::to_string(fe_degree) +
                 "):");

  auto laplace_matrix_renumbered =
    TensorProductMatrixCreator::create_1d_cell_laplace_matrix(
      fe_q, 1., {true, true}, fe_eval.get_internal_dof_numbering());
  print_matrix(laplace_matrix_renumbered,
               "Laplace matrix for FE_Q<1>(" + std::to_string(fe_degree) +
                 ") renumbered with FEEvaluation:");
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  test<1>();
  test<2>();
  test<3>();
}
