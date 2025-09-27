// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check TrilinosWrappers::SolverBelos for GMRES.

#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  using Number     = double;
  using VectorType = Vector<Number>;

  const unsigned int dim       = 2;
  const unsigned int fe_degree = 1;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  FE_Q<dim>      fe(fe_degree);
  QGauss<dim>    quad(fe_degree + 1);
  MappingQ1<dim> mapping;

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> affine_constraints;
  DoFTools::make_zero_boundary_constraints(dof_handler, affine_constraints);
  affine_constraints.close();

  TrilinosWrappers::SparsityPattern dsp(dof_handler.locally_owned_dofs(),
                                        dof_handler.get_mpi_communicator());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, affine_constraints);
  dsp.compress();

  TrilinosWrappers::SparseMatrix system_matrix;
  system_matrix.reinit(dsp);

  MatrixCreator::
    create_laplace_matrix<dim, dim, TrilinosWrappers::SparseMatrix>(
      mapping, dof_handler, quad, system_matrix, nullptr, affine_constraints);

  TrilinosWrappers::PreconditionILU ilu;
  ilu.initialize(system_matrix);

  VectorType x(dof_handler.n_dofs());
  VectorType r(dof_handler.n_dofs());

  VectorTools::create_right_hand_side(mapping,
                                      dof_handler,
                                      quad,
                                      Functions::ConstantFunction<dim, Number>(
                                        1.0),
                                      r,
                                      affine_constraints);

  Teuchos::RCP<Teuchos::ParameterList> belos_parameters =
    Teuchos::rcp(new Teuchos::ParameterList);

  belos_parameters->set("Num Blocks", 20);
  belos_parameters->set("Block Size", 10);
  belos_parameters->set("Verbosity", 0);

  if (true)
    {
      x = 0.0;

      using MV = Epetra_MultiVector;
      using OP = Epetra_Operator;

      Teuchos::RCP<Epetra_CrsMatrix> A =
        Teuchos::rcp(const_cast<Epetra_CrsMatrix *>(
                       &system_matrix.trilinos_matrix()),
                     false);
      Teuchos::RCP<Epetra_MultiVector> B, X;

      LinearAlgebra::EpetraWrappers::Vector x_(
        dof_handler.locally_owned_dofs(), dof_handler.get_mpi_communicator());
      LinearAlgebra::ReadWriteVector<Number> x_temp(
        dof_handler.locally_owned_dofs());
      x_temp.import_elements(x, VectorOperation::insert);
      x_.import_elements(x_temp, VectorOperation::insert);

      LinearAlgebra::EpetraWrappers::Vector r_(
        dof_handler.locally_owned_dofs(), dof_handler.get_mpi_communicator());
      LinearAlgebra::ReadWriteVector<Number> r_temp(
        dof_handler.locally_owned_dofs());
      r_temp.import_elements(r, VectorOperation::insert);
      r_.import_elements(r_temp, VectorOperation::insert);

      X = Teuchos::rcp<Epetra_MultiVector>(&x_.trilinos_vector(), false);
      B = Teuchos::rcp<Epetra_MultiVector>(&r_.trilinos_vector(), false);

      Belos::LinearProblem<double, MV, OP> problem(A, X, B);
      bool                                 set = problem.setProblem();

      AssertThrow(set, ExcInternalError());

      Teuchos::RCP<Belos::SolverManager<double, MV, OP>> newSolver =
        Teuchos::rcp(new Belos::BlockGmresSolMgr<double, MV, OP>(
          Teuchos::rcp(&problem, false), belos_parameters));
      Belos::ReturnType flag = newSolver->solve();

      AssertThrow(flag == Belos::ReturnType::Converged, ExcInternalError());

      deallog << x_.l2_norm() << std::endl;
    }

  if (true)
    {
      x = 0.0;

      SolverControl solver_control;
      typename TrilinosWrappers::SolverBelos<VectorType>::AdditionalData
        additional_data;

      additional_data.solver_name =
        TrilinosWrappers::SolverBelos<VectorType>::SolverName::gmres;
      additional_data.right_preconditioning = false;

      TrilinosWrappers::SolverBelos<VectorType> solver(solver_control,
                                                       additional_data,
                                                       belos_parameters);
      solver.solve(system_matrix, x, r, ilu);

      deallog << x.l2_norm() << std::endl;
    }
}
