/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2015 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------

 * This file tests the non-symmetric interface to PARPACK for an
 advection-diffusion
 * operator with PETSc mpi vectors.
 *
 * We test that the computed vectors are eigenvectors and mass-normal, i.e.
 * a) (A*x_i-\lambda*B*x_i).L2() == 0
 * b) x_i*B*x_i = 1
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
#include <deal.II/lac/parpack_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

const unsigned int dim = 2; // run in 2d to save time

const double eps = 1e-10;

template <int dim>
std::vector<dealii::IndexSet>
locally_owned_dofs_per_subdomain(const DoFHandler<dim> &dof_handler)
{
  std::vector<dealii::types::subdomain_id> subdomain_association(
    dof_handler.n_dofs());
  dealii::DoFTools::get_subdomain_association(dof_handler,
                                              subdomain_association);

  const unsigned int n_subdomains =
    1 +
    (*max_element(subdomain_association.begin(), subdomain_association.end()));

  std::vector<dealii::IndexSet> index_sets(
    n_subdomains, dealii::IndexSet(dof_handler.n_dofs()));

  // loop over subdomain_association and populate IndexSet when a
  // change in subdomain ID is found
  dealii::types::global_dof_index i_min          = 0;
  dealii::types::global_dof_index this_subdomain = subdomain_association[0];

  for (dealii::types::global_dof_index index = 1;
       index < subdomain_association.size();
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

  for (unsigned int i = 0; i < n_subdomains; ++i)
    index_sets[i].compress();

  return index_sets;
} // locally_owned_dofs_per_subdomain

class PETScInverse
{
public:
  PETScInverse(const dealii::PETScWrappers::MatrixBase &A,
               dealii::SolverControl                   &cn,
               const MPI_Comm mpi_communicator = PETSC_COMM_SELF)
    : solver(cn)
    , matrix(A)
    , preconditioner(matrix)
  {}

  void
  vmult(dealii::PETScWrappers::MPI::Vector       &dst,
        const dealii::PETScWrappers::MPI::Vector &src) const
  {
    solver.solve(matrix, dst, src, preconditioner);
  }


private:
  mutable dealii::PETScWrappers::SolverGMRES solver;
  const dealii::PETScWrappers::MatrixBase   &matrix;
  PETScWrappers::PreconditionBlockJacobi     preconditioner;
};

void
test()
{
  const unsigned int global_mesh_refinement_steps = 5;
  const unsigned int number_of_eigenvalues        = 4;

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
  std::vector<dealii::PETScWrappers::MPI::Vector> arpack_vectors;
  std::vector<std::complex<PetscScalar>>          eigenvalues;
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
    locally_owned_dofs_per_subdomain(dof_handler);
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

  // Initialise the stiffness and mass matrices
  stiffness_matrix.reinit(locally_owned_dofs,
                          locally_owned_dofs,
                          csp,
                          mpi_communicator);

  mass_matrix.reinit(locally_owned_dofs,
                     locally_owned_dofs,
                     csp,
                     mpi_communicator);

  eigenvalues.resize(number_of_eigenvalues);

  arpack_vectors.resize(number_of_eigenvalues + 1);
  for (unsigned int i = 0; i < arpack_vectors.size(); ++i)
    arpack_vectors[i].reinit(locally_owned_dofs,
                             mpi_communicator); // without ghost dofs

  eigenfunctions.resize(2 * number_of_eigenvalues);
  for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
    eigenfunctions[i].reinit(locally_owned_dofs,
                             mpi_communicator); // without ghost dofs


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
          {
            const Point<dim> cur_point = fe_values.quadrature_point(q_point);
            Tensor<1, dim>   advection;
            advection[0] = 10.;
            advection[1] = 10. * cur_point[0];
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  cell_stiffness_matrix(i, j) +=
                    (fe_values.shape_grad(i, q_point) *
                       fe_values.shape_grad(j, q_point) +
                     (advection * fe_values.shape_grad(i, q_point)) *
                       fe_values.shape_value(j, q_point)) *
                    fe_values.JxW(q_point);

                  cell_mass_matrix(i, j) +=
                    (fe_values.shape_value(i, q_point) *
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
      }

  stiffness_matrix.compress(dealii::VectorOperation::add);
  mass_matrix.compress(dealii::VectorOperation::add);

  // test Arpack
  {
    std::vector<std::complex<double>> lambda(eigenfunctions.size());

    for (unsigned int i = 0; i < eigenvalues.size(); ++i)
      eigenfunctions[i] = PetscScalar();

    dealii::SolverControl solver_control(dof_handler.n_dofs(),
                                         1e-10,
                                         /*log_history*/ false,
                                         /*log_results*/ false);
    PETScInverse inverse(stiffness_matrix, solver_control, mpi_communicator);
    const unsigned int num_arnoldi_vectors = 2 * eigenvalues.size() + 2;

    dealii::PArpackSolver<dealii::PETScWrappers::MPI::Vector>::AdditionalData
      additional_data(num_arnoldi_vectors,
                      dealii::PArpackSolver<
                        dealii::PETScWrappers::MPI::Vector>::largest_real_part,
                      false);

    dealii::PArpackSolver<dealii::PETScWrappers::MPI::Vector> eigensolver(
      solver_control, mpi_communicator, additional_data);
    eigensolver.reinit(locally_owned_dofs);
    arpack_vectors[0] = 1.;
    eigensolver.set_initial_vector(arpack_vectors[0]);
    eigensolver.solve(stiffness_matrix,
                      mass_matrix,
                      inverse,
                      eigenvalues,
                      arpack_vectors,
                      eigenvalues.size());

    // extract real and complex components of eigenvectors
    for (unsigned int i = 0; i < eigenvalues.size(); ++i)
      {
        eigenfunctions[i] = arpack_vectors[i];
        if (eigenvalues[i].imag() != 0.)
          {
            eigenfunctions[i + eigenvalues.size()] = arpack_vectors[i + 1];
            if (i + 1 < eigenvalues.size())
              {
                eigenfunctions[i + 1] = arpack_vectors[i];
                eigenfunctions[i + 1 + eigenvalues.size()] =
                  arpack_vectors[i + 1];
                eigenfunctions[i + 1 + eigenvalues.size()] *= -1;
                ++i;
              }
          }
      }

    for (unsigned int i = 0; i < eigenvalues.size(); ++i)
      dealii::deallog << eigenvalues[i] << std::endl;

    // make sure that we have eigenvectors and they are mass-normal:
    // a) (A*x_i-\lambda*B*x_i).L2() == 0
    // b) x_i*B*x_i=1
    {
      const double               precision = 1e-7;
      PETScWrappers::MPI::Vector Ax(eigenfunctions[0]), Bx(eigenfunctions[0]);
      PETScWrappers::MPI::Vector Ay(eigenfunctions[0]), By(eigenfunctions[0]);
      for (unsigned int i = 0; i < eigenvalues.size(); ++i)
        {
          stiffness_matrix.vmult(Ax, eigenfunctions[i]);
          stiffness_matrix.vmult(Ay, eigenfunctions[i + eigenvalues.size()]);
          mass_matrix.vmult(Bx, eigenfunctions[i]);
          mass_matrix.vmult(By, eigenfunctions[i + eigenvalues.size()]);

          Ax.add(-1.0 * std::real(eigenvalues[i]), Bx);
          Ax.add(std::imag(eigenvalues[i]), By);
          Ay.add(-1.0 * std::real(eigenvalues[i]), By);
          Ay.add(-1.0 * std::imag(eigenvalues[i]), Bx);
          PETScWrappers::MPI::Vector tmpx(Ax), tmpy(Ay);
          tmpx.scale(Ax);
          tmpy.scale(Ay);
          tmpx += tmpy;
          if (std::sqrt(tmpx.l1_norm()) > precision)
            deallog << "Returned vector " << i << " is not an eigenvector!"
                    << " L2 norm of the residual is "
                    << std::sqrt(tmpx.l1_norm()) << std::endl;

          const double tmp =
            std::abs(eigenfunctions[i] * Bx +
                     eigenfunctions[i + eigenvalues.size()] * By - 1.) +
            std::abs(eigenfunctions[i + eigenvalues.size()] * Bx -
                     eigenfunctions[i] * By);
          if (tmp > precision)
            deallog << "Eigenvector " << i << " is not normal! failing norm is "
                    << tmp << std::endl;
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
  deallog << std::setprecision(7);

  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);
      {
        test();
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
