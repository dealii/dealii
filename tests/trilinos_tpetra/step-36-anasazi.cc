// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Solve for smallest eigenvalues of Laplace operator in Step-36 tutorial

#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/anasazi_solver.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

namespace Step36
{
  using namespace dealii;


  template <int dim>
  class EigenvalueProblem
  {
  public:
    EigenvalueProblem();
    void
    run();

  private:
    void
    make_grid_and_dofs();
    void
    assemble_system();
    void
    solve();
    void
    output_results() const;

    Triangulation<dim> triangulation;
    const FE_Q<dim>    fe;
    DoFHandler<dim>    dof_handler;

    LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Host>
      stiffness_matrix, mass_matrix;
    std::vector<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Host>>
                        eigenfunctions;
    std::vector<double> eigenvalues;

    AffineConstraints<double> constraints;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
  };



  template <int dim>
  EigenvalueProblem<dim>::EigenvalueProblem()
    : fe(1)
    , dof_handler(triangulation)
  {}



  template <int dim>
  void
  EigenvalueProblem<dim>::make_grid_and_dofs()
  {
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(5);
    dof_handler.distribute_dofs(fe);

    DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
    constraints.close();

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    DynamicSparsityPattern dsp(locally_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               MPI_COMM_WORLD,
                                               locally_relevant_dofs);


    stiffness_matrix.reinit(locally_owned_dofs, dsp, MPI_COMM_WORLD);

    mass_matrix.reinit(locally_owned_dofs, dsp, MPI_COMM_WORLD);

    eigenfunctions.resize(5);
    for (auto &eigenfunction : eigenfunctions)
      eigenfunction.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    eigenvalues.resize(eigenfunctions.size());
  }



  template <int dim>
  void
  EigenvalueProblem<dim>::assemble_system()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
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
  }



  template <int dim>
  void
  EigenvalueProblem<dim>::solve()
  {
    SolverControl solver_control(dof_handler.n_dofs(), 1e-9);
    Teuchos::RCP<Teuchos::ParameterList> parameters =
      Teuchos::rcp(new Teuchos::ParameterList);
    parameters->set("Which", "SR");
    parameters->set("Maximum Restarts", 2000);
    parameters->set("Convergence Tolerance", 1e-9);
    parameters->set("Relative Convergence Tolerance", true);
    parameters->set("Block Size", 6);

    LinearAlgebra::TpetraWrappers::AnasaziSolverBlockDavidson<double,
                                                              MemorySpace::Host>
      eigensolver(parameters);


    eigensolver.is_hermitian(true);

    eigensolver.generate_random_initial_space(stiffness_matrix);
    eigensolver.solve(stiffness_matrix,
                      mass_matrix,
                      eigenvalues,
                      eigenfunctions,
                      eigenfunctions.size());

    for (auto &eigenfunction : eigenfunctions)
      eigenfunction /= eigenfunction.linfty_norm();
  }



  template <int dim>
  void
  EigenvalueProblem<dim>::run()
  {
    make_grid_and_dofs();

    assemble_system();

    solve();

    for (unsigned int i = 0; i < eigenvalues.size(); ++i)
      deallog << eigenvalues[i] << std::endl;

    deallog << "OK" << std::endl;
  }
} // namespace Step36

int
main(int argc, char **argv)
{
  using namespace dealii;
  using namespace Step36;

  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


  AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
              ExcMessage(
                "This program can only be run in serial, use ./step-36"));

  EigenvalueProblem<2> problem;
  problem.run();


  return 0;
}
