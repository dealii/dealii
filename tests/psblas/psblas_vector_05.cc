// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
#include <deal.II/base/exception_macros.h>
#include <deal.II/base/logstream.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/psblas_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_operation.h>

#include <psb_c_dbase.h>

#include <iostream>

#include "../tests.h"

using namespace dealii;

// Test we can perform log2glb operations with PSBLAS vectors, getting the same
// result as with PETSc vectors.

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPI_Comm                         mpi_communicator = MPI_COMM_WORLD;

  initlog();

  const unsigned int dim = 2;

  // Create distributed triangulation of the unit square
  Triangulation<dim, dim> tria_base;

  GridGenerator::hyper_cube(tria_base, 0, 1);
  tria_base.refine_global(2);

  GridTools::partition_triangulation(
    Utilities::MPI::n_mpi_processes(mpi_communicator), tria_base);

  const TriangulationDescription::Description<dim, dim> description =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria_base, mpi_communicator);

  // Create a fully distributed triangulation:
  parallel::fullydistributed::Triangulation<dim, dim> triangulation(
    mpi_communicator);
  triangulation.create_triangulation(description);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();

  AffineConstraints<double> constraints;

  PSCToolkitWrappers::Vector psblas_rhs_vector(locally_owned_dofs,
                                               mpi_communicator);
  PETScWrappers::MPI::Vector petsc_test_vector;
  petsc_test_vector.reinit(locally_owned_dofs, mpi_communicator);

  QGauss<dim>   quadrature_formula(fe.degree + 1);
  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  Vector<double>                       cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          cell_rhs = 0.;

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              const double rhs_value = 1.0;

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  cell_rhs(i) += rhs_value * fe_values.shape_value(i, q_point) *
                                 fe_values.JxW(q_point);
                }
            }

          cell->get_dof_indices(local_dof_indices);

          constraints.distribute_local_to_global(cell_rhs,
                                                 local_dof_indices,
                                                 psblas_rhs_vector);
          constraints.distribute_local_to_global(cell_rhs,
                                                 local_dof_indices,
                                                 petsc_test_vector);
        }
    }
  petsc_test_vector.compress(VectorOperation::add);
  psblas_rhs_vector.compress(VectorOperation::add);

  double difference = 0.0;
  for (const types::global_dof_index idx : locally_owned_dofs)
    difference += std::fabs(petsc_test_vector(idx) - psblas_rhs_vector(idx));

  AssertThrow(Utilities::MPI::sum(difference, mpi_communicator) < 1e-15,
              ExcInternalError());
  deallog << "OK" << std::endl;

  return 0;
}
