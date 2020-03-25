/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2018 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * This file tests the PARPACK interface for a symmetric operator taken from
 step-36
 * using Trilinos mpi vectors.
 *
 * We test that the computed vectors are eigenvectors and mass-orthonormal, i.e.
 * a) (A*x_i-\lambda*B*x_i).L2() == 0
 * b) x_j*B*x_i = \delta_{i,j}
 *
 */

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
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/parpack_solver.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

// test Parpack on Step-36 with Trilinos algebra

const unsigned int dim = 2; // run in 2d to save time


const double eps = 1e-10;

template <typename DoFHandlerType>
std::vector<IndexSet>
locally_owned_dofs_per_subdomain(const DoFHandlerType &dof_handler)
{
  std::vector<types::subdomain_id> subdomain_association(dof_handler.n_dofs());
  DoFTools::get_subdomain_association(dof_handler, subdomain_association);

  const unsigned int n_subdomains =
    1 +
    (*max_element(subdomain_association.begin(), subdomain_association.end()));

  std::vector<IndexSet> index_sets(n_subdomains,
                                   IndexSet(dof_handler.n_dofs()));

  // loop over subdomain_association and populate IndexSet when a
  // change in subdomain ID is found
  types::global_dof_index i_min          = 0;
  types::global_dof_index this_subdomain = subdomain_association[0];

  for (types::global_dof_index index = 1; index < subdomain_association.size();
       ++index)
    {
      // found index different from the current one
      if (subdomain_association[index] != this_subdomain)
        {
          index_sets[this_subdomain].add_range(i_min, index);
          i_min          = index;
          this_subdomain = subdomain_association[index];
        }
    }

  // the very last element is of different index
  if (i_min == subdomain_association.size() - 1)
    {
      index_sets[this_subdomain].add_index(i_min);
    }

  // otherwise there are at least two different indices
  else
    {
      index_sets[this_subdomain].add_range(i_min, subdomain_association.size());
    }

  for (unsigned int i = 0; i < n_subdomains; i++)
    index_sets[i].compress();

  return index_sets;
} // locally_owned_dofs_per_subdomain


void
test()
{
  const unsigned int global_mesh_refinement_steps = 5;
  const unsigned int number_of_eigenvalues        = 5;

  MPI_Comm           mpi_communicator = MPI_COMM_WORLD;
  const unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(mpi_communicator);
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);


  Triangulation<dim>        triangulation;
  DoFHandler<dim>           dof_handler(triangulation);
  FE_Q<dim>                 fe(1);
  AffineConstraints<double> constraints;
  IndexSet                  locally_owned_dofs;
  IndexSet                  locally_relevant_dofs;

  std::vector<TrilinosWrappers::MPI::Vector> eigenfunctions;
  std::vector<double>                        eigenvalues;
  TrilinosWrappers::SparseMatrix             stiffness_matrix, mass_matrix;

  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(global_mesh_refinement_steps);

  // we do not use metis but rather partition by hand below.
  // dealii::GridTools::partition_triangulation (n_mpi_processes,
  // triangulation);
  {
    const double x0 = -1.0;
    const double x1 = 1.0;
    const double dL = (x1 - x0) / n_mpi_processes;

    Triangulation<dim>::active_cell_iterator cell =
                                               triangulation.begin_active(),
                                             endc = triangulation.end();
    for (; cell != endc; ++cell)
      {
        const Point<dim> &center = cell->center();
        const double      x      = center[0];

        const auto id = static_cast<unsigned int>((x - x0) / dL);
        cell->set_subdomain_id(id);
      }
  }


  dof_handler.distribute_dofs(fe);
  DoFRenumbering::subdomain_wise(dof_handler);
  std::vector<IndexSet> locally_owned_dofs_per_processor =
    locally_owned_dofs_per_subdomain(dof_handler);
  locally_owned_dofs = locally_owned_dofs_per_processor[this_mpi_process];
  locally_relevant_dofs.clear();
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  DynamicSparsityPattern csp(locally_relevant_dofs);
  // Fill in ignoring all cells that are not locally owned
  DoFTools::make_sparsity_pattern(dof_handler,
                                  csp,
                                  constraints,
                                  /* keep constrained dofs */ true);
  std::vector<types::global_dof_index> n_locally_owned_dofs(n_mpi_processes);
  for (unsigned int i = 0; i < n_mpi_processes; i++)
    n_locally_owned_dofs[i] = locally_owned_dofs_per_processor[i].n_elements();

  SparsityTools::distribute_sparsity_pattern(csp,
                                             n_locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);

  // Initialise the stiffness and mass matrices
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
    eigenfunctions[i].reinit(locally_owned_dofs,
                             mpi_communicator); // without ghost dofs

  eigenvalues.resize(eigenfunctions.size());


  // ready for assembly
  stiffness_matrix = 0;
  mass_matrix      = 0;

  QGauss<dim>   quadrature_formula(2);
  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
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

  stiffness_matrix.compress(VectorOperation::add);
  mass_matrix.compress(VectorOperation::add);

  // test Arpack
  {
    const double                      shift = 4.0;
    std::vector<std::complex<double>> lambda(eigenfunctions.size());

    for (unsigned int i = 0; i < eigenvalues.size(); i++)
      eigenfunctions[i] = 0.;

    static ReductionControl inner_control_c(/*maxiter*/ stiffness_matrix.m(),
                                            /*tolerance (global)*/ 0.0,
                                            /*reduce (w.r.t. initial)*/ 1.e-13);

    typedef TrilinosWrappers::MPI::Vector  VectorType;
    SolverCG<VectorType>                   solver_c(inner_control_c);
    TrilinosWrappers::PreconditionIdentity preconditioner;

    const auto shifted_matrix =
      linear_operator<VectorType>(stiffness_matrix) -
      shift * linear_operator<VectorType>(mass_matrix);

    const auto shift_and_invert =
      inverse_operator(shifted_matrix, solver_c, preconditioner);

    const unsigned int num_arnoldi_vectors = 2 * eigenvalues.size() + 2;

    PArpackSolver<TrilinosWrappers::MPI::Vector>::AdditionalData
      additional_data(
        num_arnoldi_vectors,
        PArpackSolver<TrilinosWrappers::MPI::Vector>::largest_magnitude,
        true);

    SolverControl solver_control(dof_handler.n_dofs(),
                                 1e-9,
                                 /*log_history*/ false,
                                 /*log_results*/ false);

    PArpackSolver<TrilinosWrappers::MPI::Vector> eigensolver(solver_control,
                                                             mpi_communicator,
                                                             additional_data);
    eigensolver.reinit(locally_owned_dofs);
    eigensolver.set_shift(shift);
    eigenfunctions[0] = 1.;
    eigensolver.set_initial_vector(eigenfunctions[0]);
    // avoid output of iterative solver:
    const unsigned int previous_depth = deallog.depth_file(0);
    eigensolver.solve(stiffness_matrix,
                      mass_matrix,
                      shift_and_invert,
                      lambda,
                      eigenfunctions,
                      eigenvalues.size());
    deallog.depth_file(previous_depth);

    for (unsigned int i = 0; i < lambda.size(); i++)
      eigenvalues[i] = lambda[i].real();

    for (unsigned int i = 0; i < eigenvalues.size(); i++)
      deallog << eigenvalues[i] << std::endl;

    // make sure that we have eigenvectors and they are mass-orthonormal:
    // a) (A*x_i-\lambda*B*x_i).L2() == 0
    // b) x_j*B*x_i=\delta_{ij}
    {
      const double                  precision = 1e-7;
      TrilinosWrappers::MPI::Vector Ax(eigenfunctions[0]),
        Bx(eigenfunctions[0]);
      for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
        {
          mass_matrix.vmult(Bx, eigenfunctions[i]);

          for (unsigned int j = 0; j < eigenfunctions.size(); j++)
            Assert(std::abs(eigenfunctions[j] * Bx - (i == j)) < precision,
                   ExcMessage("Eigenvectors " + Utilities::int_to_string(i) +
                              " and " + Utilities::int_to_string(j) +
                              " are not orthonormal!"));

          stiffness_matrix.vmult(Ax, eigenfunctions[i]);
          Ax.add(-1.0 * std::real(lambda[i]), Bx);
          Assert(Ax.l2_norm() < precision,
                 ExcMessage("Returned vector " + Utilities::int_to_string(i) +
                            " is not an eigenvector!"));
        }
    }
  }


  dof_handler.clear();
  deallog << "Ok" << std::endl;
}


int
main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile, /*do not print job id*/ false);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        test();
      }
    }
  catch (std::exception &exc)
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
