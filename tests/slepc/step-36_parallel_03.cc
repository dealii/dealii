// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// same as step-36_parallel_02, but solve SHEP

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

// test parallel (MPI) version of Step-36

const unsigned int dim = 2; // run in 2d to save time

const double eps = 1e-10;

void
test(std::string solver_name, std::string preconditioner_name)
{
  const unsigned int global_mesh_refinement_steps = 5;
  const unsigned int number_of_eigenvalues        = 5;

  MPI_Comm           mpi_communicator = MPI_COMM_WORLD;
  const unsigned int n_mpi_processes =
    dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);
  const unsigned int this_mpi_process =
    dealii::Utilities::MPI::this_mpi_process(mpi_communicator);


  dealii::Triangulation<dim>        triangulation;
  dealii::DoFHandler<dim>           dof_handler(triangulation);
  dealii::FE_Q<dim>                 fe(1);
  dealii::AffineConstraints<double> constraints;
  dealii::IndexSet                  locally_owned_dofs;
  dealii::IndexSet                  locally_relevant_dofs;

  std::vector<dealii::PETScWrappers::MPI::Vector> eigenfunctions;
  std::vector<PetscScalar>                        eigenvalues;
  dealii::PETScWrappers::MPI::SparseMatrix        stiffness_matrix, mass_matrix;

  dealii::GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(global_mesh_refinement_steps);

  // we do not use metis but rather partition by hand below.
  // dealii::GridTools::partition_triangulation (n_mpi_processes,
  // triangulation);
  {
    const double x0 = -1.0;
    const double x1 = 1.0;
    const double dL = (x1 - x0) / n_mpi_processes;

    dealii::Triangulation<dim>::active_cell_iterator cell = triangulation
                                                              .begin_active(),
                                                     endc = triangulation.end();
    for (; cell != endc; ++cell)
      {
        const dealii::Point<dim> &center = cell->center();
        const double              x      = center[0];

        const auto id = static_cast<unsigned int>((x - x0) / dL);
        cell->set_subdomain_id(id);
      }
  }

  dof_handler.distribute_dofs(fe);
  dealii::DoFRenumbering::subdomain_wise(dof_handler);
  std::vector<dealii::IndexSet> locally_owned_dofs_per_processor =
    DoFTools::locally_owned_dofs_per_subdomain(dof_handler);
  locally_owned_dofs = locally_owned_dofs_per_processor[this_mpi_process];
  locally_relevant_dofs =
    dealii::DoFTools::extract_locally_relevant_dofs(dof_handler);

  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  dealii::VectorTools::interpolate_boundary_values(
    dof_handler, 0, dealii::Functions::ZeroFunction<dim>(), constraints);
  constraints.close();

  dealii::DynamicSparsityPattern csp(locally_relevant_dofs);
  // Fill in ignoring all cells that are not locally owned
  dealii::DoFTools::make_sparsity_pattern(dof_handler,
                                          csp,
                                          constraints,
                                          /* keep constrained dofs */ true);
  std::vector<dealii::types::global_dof_index> n_locally_owned_dofs(
    n_mpi_processes);
  for (unsigned int i = 0; i < n_mpi_processes; ++i)
    n_locally_owned_dofs[i] = locally_owned_dofs_per_processor[i].n_elements();

  dealii::SparsityTools::distribute_sparsity_pattern(csp,
                                                     n_locally_owned_dofs,
                                                     mpi_communicator,
                                                     locally_relevant_dofs);

  // initialize the stiffness and mass matrices
  stiffness_matrix.reinit(locally_owned_dofs,
                          locally_owned_dofs,
                          csp,
                          mpi_communicator);

  mass_matrix.reinit(locally_owned_dofs,
                     locally_owned_dofs,
                     csp,
                     mpi_communicator);

  eigenfunctions.resize(5);
  for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
    {
      eigenfunctions[i].reinit(locally_owned_dofs,
                               mpi_communicator); // without ghost dofs
      for (unsigned int j = 0; j < locally_owned_dofs.n_elements(); ++j)
        eigenfunctions[i][locally_owned_dofs.nth_index_in_set(j)] =
          random_value<double>();

      eigenfunctions[i].compress(dealii::VectorOperation::insert);
    }

  eigenvalues.resize(eigenfunctions.size());


  // ready for assembly
  stiffness_matrix = 0;
  mass_matrix      = 0;

  dealii::QGauss<dim>   quadrature_formula(2);
  dealii::FEValues<dim> fe_values(fe,
                                  quadrature_formula,
                                  dealii::update_values |
                                    dealii::update_gradients |
                                    dealii::update_quadrature_points |
                                    dealii::update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  dealii::FullMatrix<double> cell_stiffness_matrix(dofs_per_cell,
                                                   dofs_per_cell);
  dealii::FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

  typename dealii::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell != endc; ++cell)
    if (cell->subdomain_id() == this_mpi_process)
      {
        fe_values.reinit(cell);
        cell_stiffness_matrix = 0;
        cell_mass_matrix      = 0;

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                cell_stiffness_matrix(i, j) +=
                  (fe_values.shape_grad(i, q_point) *
                   fe_values.shape_grad(j, q_point)) *
                  fe_values.JxW(q_point);

                cell_mass_matrix(i, j) += (fe_values.shape_value(i, q_point) *
                                           fe_values.shape_value(j, q_point)) *
                                          fe_values.JxW(q_point);
              }

        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(cell_stiffness_matrix,
                                               local_dof_indices,
                                               stiffness_matrix);
        constraints.distribute_local_to_global(cell_mass_matrix,
                                               local_dof_indices,
                                               mass_matrix);
      }

  stiffness_matrix.compress(dealii::VectorOperation::add);
  mass_matrix.compress(dealii::VectorOperation::add);

  // test SLEPc by
  {
    PETScWrappers::PreconditionBase *preconditioner = nullptr;

    dealii::deallog << preconditioner_name << std::endl;
    if (preconditioner_name == "Jacobi")
      {
        preconditioner =
          new PETScWrappers::PreconditionJacobi(mpi_communicator);
      }
    else if (preconditioner_name == "Boomer")
      {
        PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
        data.symmetric_operator = true;

        preconditioner =
          new PETScWrappers::PreconditionBoomerAMG(mpi_communicator, data);
      }
    else if (preconditioner_name == "BlockJacobi")
      {
        preconditioner =
          new PETScWrappers::PreconditionBlockJacobi(mpi_communicator);
      }
    else
      {
        AssertThrow(false, ExcMessage("Unsupported preconditioner"));
      }

    dealii::SolverControl   linear_solver_control(dof_handler.n_dofs(),
                                                1e-15,
                                                /*log_history*/ false,
                                                /*log_results*/ false);
    PETScWrappers::SolverCG linear_solver(linear_solver_control);
    linear_solver.initialize(*preconditioner);

    dealii::SolverControl solver_control(100,
                                         1e-12,
                                         /*log_history*/ false,
                                         /*log_results*/ false);

    dealii::SLEPcWrappers::SolverBase *eigensolver;

    dealii::deallog << solver_name << std::endl;
    // Get a handle on the wanted eigenspectrum solver
    if (solver_name == "KrylovSchur")
      {
        eigensolver =
          new dealii::SLEPcWrappers::SolverKrylovSchur(solver_control,
                                                       mpi_communicator);
      }

    else if (solver_name == "GeneralizedDavidson")
      {
        eigensolver = new dealii::SLEPcWrappers::SolverGeneralizedDavidson(
          solver_control, mpi_communicator);
      }
    else if (solver_name == "JacobiDavidson")
      {
        eigensolver =
          new dealii::SLEPcWrappers::SolverJacobiDavidson(solver_control,
                                                          mpi_communicator);
      }
    else if (solver_name == "Lanczos")
      {
        eigensolver =
          new dealii::SLEPcWrappers::SolverLanczos(solver_control,
                                                   mpi_communicator);
      }
    else
      {
        AssertThrow(false, ExcMessage("not supported eigensolver"));

        // Make compiler happy and not complaining about non
        // uninitialized variables
        eigensolver =
          new dealii::SLEPcWrappers::SolverKrylovSchur(solver_control,
                                                       mpi_communicator);
      }

    // Set the initial vector. This is optional, if not done the initial vector
    // is set to random values
    eigensolver->set_initial_space(eigenfunctions);

    eigensolver->set_which_eigenpairs(EPS_LARGEST_REAL);
    eigensolver->set_problem_type(EPS_HEP);

    eigensolver->solve(stiffness_matrix,
                       eigenvalues,
                       eigenfunctions,
                       eigenfunctions.size());

    // TODO make this robust on different platforms. Seems related to GHEP
    // as solve_04 works ok.
    // dealii::deallog << "outer iterations: "<< solver_control.last_step
    // ()<<std::endl; dealii::deallog << "last inner iterations:
    // "<<linear_solver_control.last_step()<<std::endl;
    for (unsigned int i = 0; i < eigenvalues.size(); ++i)
      dealii::deallog << eigenvalues[i] << std::endl;

    delete preconditioner;
    delete eigensolver;

    // make sure that we have eigenvectors and they are mass-orthonormal:
    // a) (A*x_i-\lambda*x_i).L2() == 0
    // b) x_j*x_i=\delta_{ij}
    {
      const double               precision = 1e-5;
      PETScWrappers::MPI::Vector Ax(eigenfunctions[0]);
      for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
        {
          for (unsigned int j = 0; j < eigenfunctions.size(); ++j)
            Assert(std::abs(eigenfunctions[j] * eigenfunctions[i] - (i == j)) <
                     precision,
                   ExcMessage("Eigenvectors " + Utilities::int_to_string(i) +
                              " and " + Utilities::int_to_string(j) +
                              " are not orthonormal!"));

          stiffness_matrix.vmult(Ax, eigenfunctions[i]);
          Ax.add(-1.0 * eigenvalues[i], eigenfunctions[i]);
          Assert(Ax.l2_norm() < precision,
                 ExcMessage(Utilities::to_string(Ax.l2_norm())));
        }
    }
  }


  dof_handler.clear();
  dealii::deallog << "Ok" << std::endl;
}


int
main(int argc, char **argv)
{
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile, /*do not print job id*/ false);
  dealii::deallog.depth_console(0);

  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);
      {
        test("KrylovSchur", "Jacobi");
        test("KrylovSchur", "BlockJacobi");
        test("KrylovSchur", "Boomer");
        //        test ("GeneralizedDavidson");
        //        test ("JacobiDavidson");
      }
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
