// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test the mumps sparse direct solver on a mass matrix
// - check several ways to solve, also using block low rank preconditioning

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "../tests.h"


void
solve_and_check(const SparseMatrix<double> &M,
                const Vector<double>       &rhs,
                const Vector<double>       &solution)
{
  {
    int                              *icntl = nullptr;
    SparseDirectMUMPS::AdditionalData data;
    data.output_details = false;
    SparseDirectMUMPS solver(data);
    icntl = solver.get_icntl();
    std::cout << "Elementi di icntl " << icntl[7] << " " << icntl[8] << " "
              << icntl[9] << std::endl;
    solver.initialize(M);
    Vector<double> dst(rhs.size());
    solver.vmult(dst, rhs);
    dst -= solution;
    Assert(dst.l2_norm() < 1e-9, ExcInternalError());
  }
}



template <int dim>
void
test()
{
  deallog << dim << 'd' << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1);
  tria.refine_global(1);

  // destroy the uniformity of the matrix by
  // refining one cell
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(8 - 2 * dim);

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

  SparseMatrix<double> B;
  B.reinit(sparsity_pattern);

  QGauss<dim> qr(2);
  MatrixTools::create_mass_matrix(dof_handler, qr, B);

  // compute a decomposition of the matrix
  SparseDirectMUMPS::AdditionalData data;
  data.blr_factorization = true;
  SparseDirectMUMPS::AdditionalData::BlockLowRank blr;
  blr.blr_ucfs          = false;
  blr.lowrank_threshold = 1e-8;
  data.output_details   = false;
  data.error_statistics = false;
  SparseDirectMUMPS Binv(data);

  Binv.initialize(B);

  // for a number of different solution
  // vectors, make up a matching rhs vector
  // and check what the solver finds
  for (unsigned int i = 0; i < 3; ++i)
    {
      Vector<double> solution(dof_handler.n_dofs());
      Vector<double> x(dof_handler.n_dofs());
      Vector<double> b(dof_handler.n_dofs());

      for (unsigned int j = 0; j < dof_handler.n_dofs(); ++j)
        solution(j) = j + j * (i + 1) * (i + 1);

      B.vmult(b, solution);

      SolverControl            solver_control(1000, 1e-8);
      SolverCG<Vector<double>> solver(solver_control);

      solver.solve(B, x, b, Binv);

      deallog << "Number of preconditioned CG iterations = "
              << solver_control.last_step() << std::endl;

      x -= solution;
      deallog << "relative norm distance = " << x.l2_norm() / solution.l2_norm()
              << std::endl;
      deallog << "absolute norms = " << x.l2_norm() << ' ' << solution.l2_norm()
              << std::endl;
      Assert(x.l2_norm() / solution.l2_norm() < 1e-8, ExcInternalError());

      solve_and_check(B, b, solution);
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  initlog();

  test<1>();
  test<2>();
  test<3>();
}
