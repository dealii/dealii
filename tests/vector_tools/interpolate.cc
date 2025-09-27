// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// The test initializes a non-zero solution over a hypercube using the
// finite element space FESystem(FE_Q(2), dim). The solution is then
// interpolated to the finite element space FESystem(FE_Q(3), dim) using
// the parallel VectorTools::interpolate() function.
//
// Since the solution is continuous across the elements, the solution
// values at the new support points of FESystem(FE_Q(3), dim) should be
// equal to the old solution defined in FESystem(FE_Q(2), dim) at those
// same support points.

#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
run(const dealii::Triangulation<dim> &triangulation)
{
  using Vector = dealii::LinearAlgebra::distributed::Vector<double>;

  // Create low order system.
  dealii::DoFHandler<dim> dof_handler_coarse(triangulation);

  const unsigned int          degree_coarse = 2;
  const dealii::FE_Q<dim>     fe_q_coarse(degree_coarse);
  const dealii::FESystem<dim> fe_system_coarse(fe_q_coarse, dim);
  dof_handler_coarse.distribute_dofs(fe_system_coarse);

  Vector                 solution_coarse;
  const dealii::IndexSet locally_owned_dofs_coarse =
    dof_handler_coarse.locally_owned_dofs();
  dealii::IndexSet ghost_dofs_coarse =
    dealii::DoFTools::extract_locally_relevant_dofs(dof_handler_coarse);
  dealii::IndexSet locally_relevant_dofs_coarse = ghost_dofs_coarse;
  ghost_dofs_coarse.subtract_set(locally_owned_dofs_coarse);

  if (dim == 1)
    {
      solution_coarse.reinit(locally_owned_dofs_coarse, MPI_COMM_WORLD);
    }
  else
    {
      solution_coarse.reinit(locally_owned_dofs_coarse,
                             ghost_dofs_coarse,
                             MPI_COMM_WORLD);
    }

  // Initialize dummy coarse solution.
  for (unsigned int idof = 0; idof < dof_handler_coarse.n_dofs(); ++idof)
    {
      if (locally_owned_dofs_coarse.is_element(idof))
        {
          solution_coarse[idof] = idof;
        }
    }
  solution_coarse.update_ghost_values();

  const unsigned int degree_fine = 3;

  dealii::DoFHandler<dim> dof_handler_fine(triangulation);

  const dealii::FE_Q<dim>     fe_q_fine(degree_fine);
  const dealii::FESystem<dim> fe_system_fine(fe_q_fine, dim);
  dof_handler_fine.distribute_dofs(fe_system_fine);

  Vector                 solution_fine;
  const dealii::IndexSet locally_owned_dofs_fine =
    dof_handler_fine.locally_owned_dofs();
  dealii::IndexSet ghost_dofs_fine =
    dealii::DoFTools::extract_locally_relevant_dofs(dof_handler_fine);
  dealii::IndexSet locally_relevant_dofs_fine = ghost_dofs_fine;
  ghost_dofs_fine.subtract_set(locally_owned_dofs_fine);

  if (dim == 1)
    {
      solution_fine.reinit(locally_owned_dofs_fine, MPI_COMM_WORLD);
    }
  else
    {
      solution_fine.reinit(locally_owned_dofs_fine,
                           ghost_dofs_fine,
                           MPI_COMM_WORLD);
    }

  // Interpolate the solution.
  dealii::FullMatrix<double> interpolation_matrix(
    fe_system_fine.dofs_per_cell, fe_system_coarse.dofs_per_cell);
  fe_system_fine.get_interpolation_matrix(fe_system_coarse,
                                          interpolation_matrix);
  dealii::VectorTools::interpolate(dof_handler_coarse,
                                   dof_handler_fine,
                                   interpolation_matrix,
                                   solution_coarse,
                                   solution_fine);
  solution_fine.update_ghost_values();

  const unsigned int n_dofs_coarse = fe_system_coarse.dofs_per_cell;
  const unsigned int n_dofs_fine   = fe_system_fine.dofs_per_cell;
  const std::vector<dealii::Point<dim>> &points =
    fe_system_fine.get_unit_support_points();

  // Check that interpolated solution matches at quadrature points.
  std::vector<dealii::types::global_dof_index> dof_indices_coarse(
    n_dofs_coarse);
  std::vector<dealii::types::global_dof_index> dof_indices_fine(n_dofs_fine);

  auto cell_coarse = dof_handler_coarse.begin_active();
  auto cell_fine   = dof_handler_fine.begin_active();
  auto endcell     = dof_handler_fine.end();

  bool matching = true;
  for (; cell_fine != endcell; ++cell_fine, ++cell_coarse)
    {
      if (!cell_fine->is_locally_owned())
        continue;

      cell_fine->get_dof_indices(dof_indices_fine);
      cell_coarse->get_dof_indices(dof_indices_coarse);

      for (unsigned int iquad = 0; iquad < points.size(); ++iquad)
        {
          dealii::Tensor<1, dim, double> local_solution_coarse;
          dealii::Tensor<1, dim, double> local_solution_fine;

          for (unsigned int idof = 0; idof < n_dofs_coarse; ++idof)
            {
              const unsigned int axis =
                fe_system_coarse.system_to_component_index(idof).first;
              local_solution_coarse[axis] +=
                solution_coarse[dof_indices_coarse[idof]] *
                fe_system_coarse.shape_value(idof, points[iquad]);
            }

          for (unsigned int idof = 0; idof < n_dofs_fine; ++idof)
            {
              const unsigned int axis =
                fe_system_fine.system_to_component_index(idof).first;
              local_solution_fine[axis] +=
                solution_fine[dof_indices_fine[idof]] *
                fe_system_fine.shape_value(idof, points[iquad]);
            }

          dealii::Tensor<1, dim, double> diff = local_solution_fine;
          diff -= local_solution_coarse;

          const double diff_norm = diff.norm();

          std::cout << "Coarse solution = " << local_solution_coarse
                    << ", fine solution = " << local_solution_fine
                    << ", difference = " << diff_norm << std::endl;

          if (diff_norm > 1e-12)
            {
              matching = false;
            }
        }
    }

  if (matching)
    {
      dealii::deallog << "OK" << std::endl;
    }
  else
    {
      dealii::deallog << "Non matching" << std::endl;
    }
}

int
main(int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  const unsigned int n_1d_cells = 2;

  {
    // 1D
    dealii::Triangulation<1> triangulation;
    dealii::GridGenerator::subdivided_hyper_cube(triangulation, n_1d_cells);

    dealii::deallog << "1D" << std::endl;
    run<1>(triangulation);
  }

  // 2D.
  {
    dealii::parallel::distributed::Triangulation<2> triangulation(
      MPI_COMM_WORLD);
    dealii::GridGenerator::subdivided_hyper_cube(triangulation, n_1d_cells);

    dealii::deallog << "2D" << std::endl;
    run<2>(triangulation);
  }

  // 3D.
  {
    dealii::parallel::distributed::Triangulation<3> triangulation(
      MPI_COMM_WORLD);
    dealii::GridGenerator::subdivided_hyper_cube(triangulation, n_1d_cells);

    dealii::deallog << "3D" << std::endl;
    run<3>(triangulation);
  }
}
