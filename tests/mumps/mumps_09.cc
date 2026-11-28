// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test MUMPS initialize_from_individual_blocks() for the Stokes
// (saddle-point) system, based on the Kovasznay flow problem from step-55.
// Instead of passing a monolithic PETScWrappers::MPI::BlockSparseMatrix to
// initialize(), this test passes pointers to the individual
// PETScWrappers::MPI::SparseMatrix sub-blocks into
// initialize_from_individual_blocks(). The (1,1) pressure-pressure block is
// passed as nullptr to exercise the zero-block code path.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include "../tests/tests.h"

// Same as mumps_08.cc, but with initialize_from_individual_blocks()
// function. Blocks are given separately to MUMPS.

namespace StokesTest
{
  using namespace dealii;

  // Kovasznay exact solution (2D, dim+1 components: u_x, u_y, p)
  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution()
      : Function<dim>(dim + 1)
    {}

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;
  };

  template <int dim>
  void
  ExactSolution<dim>::vector_value(const Point<dim> &p,
                                   Vector<double>   &values) const
  {
    const double R_x = p[0];
    const double R_y = p[1];

    const double pi  = numbers::PI;
    const double pi2 = pi * pi;

    values[0] =
      -exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi) + 1;
    values[1] = (1.0L / 2.0L) * (-sqrt(25.0 + 4 * pi2) + 5.0) *
                exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) /
                pi;
    values[2] =
      -1.0L / 2.0L * exp(R_x * (-2 * sqrt(25.0 + 4 * pi2) + 10.0)) -
      2.0 *
        (-6538034.74494422 +
         0.0134758939981709 * exp(4 * sqrt(25.0 + 4 * pi2))) /
        (-80.0 * exp(3 * sqrt(25.0 + 4 * pi2)) +
         16.0 * sqrt(25.0 + 4 * pi2) * exp(3 * sqrt(25.0 + 4 * pi2))) -
      1634508.68623606 * exp(-3.0 * sqrt(25.0 + 4 * pi2)) /
        (-10.0 + 2.0 * sqrt(25.0 + 4 * pi2)) +
      (-0.00673794699908547 * exp(sqrt(25.0 + 4 * pi2)) +
       3269017.37247211 * exp(-3 * sqrt(25.0 + 4 * pi2))) /
        (-8 * sqrt(25.0 + 4 * pi2) + 40.0) +
      0.00336897349954273 * exp(1.0 * sqrt(25.0 + 4 * pi2)) /
        (-10.0 + 2.0 * sqrt(25.0 + 4 * pi2));
  }


  // Kovasznay right-hand side (body force for Stokes with viscosity=0.1)
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>(dim + 1)
    {}

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;
  };

  template <int dim>
  void
  RightHandSide<dim>::vector_value(const Point<dim> &p,
                                   Vector<double>   &values) const
  {
    const double R_x = p[0];
    const double R_y = p[1];

    const double pi  = numbers::PI;
    const double pi2 = pi * pi;

    values[0] =
      -1.0L / 2.0L * (-2 * sqrt(25.0 + 4 * pi2) + 10.0) *
        exp(R_x * (-2 * sqrt(25.0 + 4 * pi2) + 10.0)) -
      0.4 * pi2 * exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi) +
      0.1 * pow(-sqrt(25.0 + 4 * pi2) + 5.0, 2) *
        exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi);
    values[1] = 0.2 * pi * (-sqrt(25.0 + 4 * pi2) + 5.0) *
                  exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) -
                0.05 * pow(-sqrt(25.0 + 4 * pi2) + 5.0, 3) *
                  exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) /
                  pi;
    values[2] = 0;
  }


  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem(unsigned int velocity_degree);
    void
    run();

  private:
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    compute_errors() const;

    const unsigned int velocity_degree;
    const double       viscosity;
    MPI_Comm           mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;

    std::vector<IndexSet> owned_partitioning;
    std::vector<IndexSet> relevant_partitioning;

    AffineConstraints<double> constraints;

    PETScWrappers::MPI::BlockSparseMatrix system_matrix;
    PETScWrappers::MPI::BlockVector       locally_relevant_solution;
    PETScWrappers::MPI::BlockVector       system_rhs;
  };


  template <int dim>
  StokesProblem<dim>::StokesProblem(unsigned int velocity_degree)
    : velocity_degree(velocity_degree)
    , viscosity(0.1)
    , mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , fe(FE_Q<dim>(velocity_degree), dim, FE_Q<dim>(velocity_degree - 1), 1)
    , dof_handler(triangulation)
  {}


  template <int dim>
  void
  StokesProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0);
    stokes_sub_blocks[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, stokes_sub_blocks);

    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, stokes_sub_blocks);

    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    deallog << "Number of dofs = " << dof_handler.n_dofs() << " (" << n_u
            << " + " << n_p << ")" << std::endl;

    owned_partitioning.resize(2);
    owned_partitioning[0] = dof_handler.locally_owned_dofs().get_view(0, n_u);
    owned_partitioning[1] =
      dof_handler.locally_owned_dofs().get_view(n_u, n_u + n_p);

    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    relevant_partitioning.resize(2);
    relevant_partitioning[0] = locally_relevant_dofs.get_view(0, n_u);
    relevant_partitioning[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    {
      constraints.reinit(dof_handler.locally_owned_dofs(),
                         locally_relevant_dofs);

      const FEValuesExtractors::Vector velocities(0);
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               ExactSolution<dim>(),
                                               constraints,
                                               fe.component_mask(velocities));
      constraints.close();
    }

    {
      system_matrix.clear();

      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
      for (unsigned int c = 0; c < dim + 1; ++c)
        for (unsigned int d = 0; d < dim + 1; ++d)
          if (c == dim && d == dim)
            coupling[c][d] = DoFTools::none;
          else if (c == dim || d == dim || c == d)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

      DoFTools::make_sparsity_pattern(
        dof_handler, coupling, dsp, constraints, false);

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        dof_handler.locally_owned_dofs(),
        mpi_communicator,
        locally_relevant_dofs);

      system_matrix.reinit(owned_partitioning, dsp, mpi_communicator);
    }

    locally_relevant_solution.reinit(owned_partitioning,
                                     relevant_partitioning,
                                     mpi_communicator);
    system_rhs.reinit(owned_partitioning, mpi_communicator);
  }


  template <int dim>
  void
  StokesProblem<dim>::assemble_system()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const QGauss<dim> quadrature_formula(velocity_degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    const RightHandSide<dim>    right_hand_side;
    std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1));

    std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
    std::vector<double>         div_phi_u(dofs_per_cell);
    std::vector<double>         phi_p(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    const FEValuesExtractors::Vector     velocities(0);
    const FEValuesExtractors::Scalar     pressure(dim);

    for (const auto &cell : dof_handler.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
      {
        cell_matrix = 0;
        cell_rhs    = 0;

        fe_values.reinit(cell);
        right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                          rhs_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                div_phi_u[k]  = fe_values[velocities].divergence(k, q);
                phi_p[k]      = fe_values[pressure].value(k, q);
              }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    cell_matrix(i, j) +=
                      (viscosity *
                         scalar_product(grad_phi_u[i], grad_phi_u[j]) -
                       div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
                      fe_values.JxW(q);
                  }

                const unsigned int component_i =
                  fe.system_to_component_index(i).first;
                cell_rhs(i) += fe_values.shape_value(i, q) *
                               rhs_values[q](component_i) * fe_values.JxW(q);
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }


  template <int dim>
  void
  StokesProblem<dim>::solve()
  {
    // Stokes system is indefinite
    SparseDirectMUMPS::AdditionalData data;
    data.output_details = false;
    data.symmetric      = false;
    data.posdef         = false;
    SparseDirectMUMPS solver(data, system_matrix.get_mpi_communicator());

    // Instead of solver.initialize(system_matrix), pass individual blocks.
    // The (1,1) pressure-pressure block is zero, so we pass nullptr.
    std::vector<const PETScWrappers::MPI::SparseMatrix *> blocks = {
      &system_matrix.block(0, 0),
      &system_matrix.block(0, 1),
      &system_matrix.block(1, 0),
      nullptr};
    solver.initialize_from_individual_blocks(blocks, 2, 2);

    PETScWrappers::MPI::BlockVector distributed_solution(owned_partitioning,
                                                         mpi_communicator);

    solver.vmult(distributed_solution, system_rhs);

    constraints.distribute(distributed_solution);

    // Subtract mean pressure to match the exact solution
    locally_relevant_solution = distributed_solution;
    const double mean_pressure =
      VectorTools::compute_mean_value(dof_handler,
                                      QGauss<dim>(velocity_degree + 2),
                                      locally_relevant_solution,
                                      dim);
    distributed_solution.block(1).add(-mean_pressure);
    locally_relevant_solution.block(1) = distributed_solution.block(1);
  }


  template <int dim>
  void
  StokesProblem<dim>::compute_errors() const
  {
    const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1);
    const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
                                                     dim + 1);

    Vector<double>    cellwise_errors(triangulation.n_active_cells());
    const QGauss<dim> quadrature(velocity_degree + 2);

    VectorTools::integrate_difference(dof_handler,
                                      locally_relevant_solution,
                                      ExactSolution<dim>(),
                                      cellwise_errors,
                                      quadrature,
                                      VectorTools::L2_norm,
                                      &velocity_mask);

    const double error_u_l2 =
      VectorTools::compute_global_error(triangulation,
                                        cellwise_errors,
                                        VectorTools::L2_norm);

    VectorTools::integrate_difference(dof_handler,
                                      locally_relevant_solution,
                                      ExactSolution<dim>(),
                                      cellwise_errors,
                                      quadrature,
                                      VectorTools::L2_norm,
                                      &pressure_mask);

    const double error_p_l2 =
      VectorTools::compute_global_error(triangulation,
                                        cellwise_errors,
                                        VectorTools::L2_norm);

    deallog << "error: u_0: " << error_u_l2 << " p_0: " << error_p_l2
            << std::endl;
  }


  template <int dim>
  void
  StokesProblem<dim>::run()
  {
    GridGenerator::hyper_cube(triangulation, -0.5, 1.5);
    triangulation.refine_global(3);

    setup_system();
    assemble_system();
    solve();
    compute_errors();
  }
} // namespace StokesTest


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log(true);

  {
    using namespace StokesTest;
    StokesProblem<2> problem(2);
    problem.run();
  }

  return 0;
}
