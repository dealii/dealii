// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Arpack by calculating eigenvalues of Laplace matrix.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/arpack_solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <complex>
#include <iostream>
#include <vector>

#include "../tests.h"


template <int dim>
void
test()
{
  const unsigned int n_eigenvalues = 3;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1);
  tria.refine_global(3);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << "Number of dofs = " << dof_handler.n_dofs() << std::endl;

  SparsityPattern sparsity_pattern;
  sparsity_pattern.reinit(dof_handler.n_dofs(),
                          dof_handler.n_dofs(),
                          dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  SparseMatrix<double> A, B;
  A.reinit(sparsity_pattern);
  B.reinit(sparsity_pattern);

  std::vector<Vector<double>> eigenvectors(
    n_eigenvalues, Vector<double>(dof_handler.n_dofs()));
  std::vector<std::complex<double>> eigenvalues(n_eigenvalues);

  QGauss<dim> qr(2);
  MatrixTools::create_laplace_matrix(dof_handler, qr, A);
  MatrixTools::create_mass_matrix(dof_handler, qr, B);

  SolverControl       solver_control(dof_handler.n_dofs(), 1e-10);
  SparseDirectUMFPACK inverse;
  inverse.initialize(A);
  const unsigned int           num_arnoldi_vectors = 2 * eigenvalues.size() + 2;
  ArpackSolver::AdditionalData additional_data(num_arnoldi_vectors,
                                               ArpackSolver::largest_magnitude,
                                               true);
  ArpackSolver                 eigensolver(solver_control, additional_data);
  eigensolver.solve(
    A, B, inverse, eigenvalues, eigenvectors, eigenvalues.size());

  {
    const double   precision = 1e-7;
    Vector<double> Ax(eigenvectors[0]), Bx(eigenvectors[0]);
    for (unsigned int i = 0; i < eigenvectors.size(); ++i)
      {
        B.vmult(Bx, eigenvectors[i]);

        for (unsigned int j = 0; j < eigenvectors.size(); ++j)
          Assert(std::abs(eigenvectors[j] * Bx - (i == j)) < precision,
                 ExcMessage("Eigenvectors " + Utilities::int_to_string(i) +
                            " and " + Utilities::int_to_string(j) +
                            " are not orthonormal!"));

        A.vmult(Ax, eigenvectors[i]);
        Ax.add(-1.0 * std::real(eigenvalues[i]), Bx);
        Assert(Ax.l2_norm() < precision,
               ExcMessage("Returned vector " + Utilities::int_to_string(i) +
                          " is not an eigenvector!"));
      }
  }
}


int
main()
{
  initlog();

  test<2>();
}
