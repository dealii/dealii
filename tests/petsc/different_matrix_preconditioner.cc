// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// If different matrices for building the preconditioner and the linear system
// to be solved are used, check whether the correct equation is solved
// https://github.com/dealii/dealii/issues/1978


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  // create stiffness and mass matrix templates

  Triangulation<2> triangulation;
  GridGenerator::hyper_cube(triangulation, 0, 1.);
  triangulation.refine_global(5);

  FE_Q<2>       fe_q(1);
  DoFHandler<2> dof_handler;
  dof_handler.initialize(triangulation, fe_q);

  QGauss<2> quadrature(2);

  SparsityPattern sparsity_pattern(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  SparseMatrix<double> A_tmp;
  A_tmp.reinit(sparsity_pattern);
  MatrixCreator::create_mass_matrix(dof_handler, quadrature, A_tmp);

  // convert to PETSc types

  PETScWrappers::SparseMatrix A;
  PETScWrappers::SparseMatrix A2;
  A.reinit(sparsity_pattern);
  A2.reinit(sparsity_pattern);

  SparseMatrix<double>::iterator it = A_tmp.begin(), endit = A_tmp.end();
  for (; it != endit; ++it)
    {
      A.set(it->row(), it->column(), it->value());
      A2.set(it->row(), it->column(), 1000. * it->value()); // scale A2 by 1000.
    }
  A.compress(VectorOperation::insert);
  A2.compress(VectorOperation::insert);

  const unsigned int n_dofs = dof_handler.n_dofs();

  PETScWrappers::MPI::Vector right_hand_side(MPI_COMM_WORLD, n_dofs, n_dofs);
  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    right_hand_side[i] = 1.;
  right_hand_side.compress(VectorOperation::insert);

  // solve
  PETScWrappers::MPI::Vector solution(MPI_COMM_WORLD, n_dofs, n_dofs);
  PETScWrappers::MPI::Vector residual(MPI_COMM_WORLD, n_dofs, n_dofs);

  {
    SolverControl              solver_control(n_dofs, 1e-6, false, false);
    PETScWrappers::SolverGMRES solver(solver_control);

    PETScWrappers::PreconditionBoomerAMG preconditioner(A);

    solver.solve(A, solution, right_hand_side, preconditioner);

    deallog << "prec A , matrix A : ";
    A.residual(residual, solution, right_hand_side);
    if (residual.l2_norm() > 1.e-6)
      deallog << residual.l2_norm() << std::endl;
    else
      deallog << "OK" << std::endl;
  }
  solution = 0.;
  {
    SolverControl              solver_control(n_dofs, 1e-6, false, false);
    PETScWrappers::SolverGMRES solver(solver_control);

    PETScWrappers::PreconditionBoomerAMG preconditioner(A2);

    solver.solve(A, solution, right_hand_side, preconditioner);

    deallog << "prec A2, matrix A, A residual : ";
    A.residual(residual, solution, right_hand_side);
    if (residual.l2_norm() > 1.e-6)
      deallog << residual.l2_norm() << std::endl;
    else
      deallog << "OK" << std::endl;

    deallog << "prec A2, matrix A, A2 residual : ";
    A2.residual(residual, solution, right_hand_side);
    if (residual.l2_norm() < 1.e-6)
      deallog << residual.l2_norm() << std::endl;
    else
      deallog << "OK" << std::endl;
  }
}
