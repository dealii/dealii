// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Check TrilinosWrappers::SolverBelos for GMRES.
// Like solver_belos_01, but with Tpetra vectors and matrices.

#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_tpetra_precondition.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <BelosTpetraAdapter.hpp>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  using Number     = double;
  using VectorType = LinearAlgebra::TpetraWrappers::Vector<Number>;
  using MatrixType = LinearAlgebra::TpetraWrappers::SparseMatrix<Number>;

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

  DynamicSparsityPattern dsp(dof_handler.locally_owned_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, affine_constraints);
  dsp.compress();

  MatrixType system_matrix;
  system_matrix.reinit(dsp);

  MatrixCreator::create_laplace_matrix<dim, dim, MatrixType>(
    mapping, dof_handler, quad, system_matrix, nullptr, affine_constraints);

  LinearAlgebra::TpetraWrappers::PreconditionILU<Number> ilu;
  ilu.initialize(system_matrix);

  VectorType x(dof_handler.locally_owned_dofs(),
               dof_handler.get_mpi_communicator());
  VectorType r(dof_handler.locally_owned_dofs(),
               dof_handler.get_mpi_communicator());

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

      using MV = LinearAlgebra::TpetraWrappers::TpetraTypes::
        MultiVectorType<Number, MemorySpace::Host>;
      using OP = LinearAlgebra::TpetraWrappers::TpetraTypes::
        LinearOperator<Number, MemorySpace::Host>;

      Teuchos::RCP<OP> A = system_matrix.trilinos_rcp();

      VectorType x_(x);
      VectorType r_(r);

      Teuchos::RCP<MV> X = Teuchos::rcp(&x_.trilinos_vector(), false);
      Teuchos::RCP<MV> B = Teuchos::rcp(&r_.trilinos_vector(), false);

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
