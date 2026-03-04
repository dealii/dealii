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

#include "deal.II/base/exceptions.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/psblas_sparse_matrix.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

// This some methods for the PSBLAS matrix class.

using namespace dealii;

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPI_Comm                         mpi_communicator = MPI_COMM_WORLD;
  MPILogInitAll                    log;

  static constexpr unsigned int dim = 2;

  // Create distributed triangulation of the unit square
  Triangulation<dim> tria_base;

  // Create a serial triangulation (here by reading an external mesh):
  GridGenerator::hyper_cube(tria_base, 0, 1);
  tria_base.refine_global(5);

  // Partition
  GridTools::partition_triangulation(
    Utilities::MPI::n_mpi_processes(mpi_communicator), tria_base);

  // Create building blocks:
  const TriangulationDescription::Description<dim> description =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria_base, mpi_communicator);

  // Create a fully distributed triangulation:
  parallel::fullydistributed::Triangulation<dim> triangulation(
    mpi_communicator);
  triangulation.create_triangulation(description);

  // Finite element and DoFHandler
  FE_Q<dim>       fe(dim);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);


  PSCToolkitWrappers::SparsityPattern sparsity_pattern(locally_owned_dofs,
                                                       mpi_communicator);
  PSCToolkitWrappers::SparseMatrix    psblas_matrix(sparsity_pattern,
                                                 mpi_communicator);

  // Assemble system
  QGauss<dim>   quadrature_formula(fe.degree + 1);
  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  AffineConstraints<double>            constraints;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell_matrix = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) += fe_values.shape_grad(i, q) *
                                     fe_values.shape_grad(j, q) *
                                     fe_values.JxW(q);

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 local_dof_indices,
                                                 psblas_matrix);
        }
    }
  psblas_matrix.compress();

  deallog << "Local range: " << psblas_matrix.local_range().first << " - "
          << psblas_matrix.local_range().second << std::endl;
  AssertThrow((psblas_matrix.in_local_range(*locally_owned_dofs.begin()) ==
               true),
              ExcMessage("First locally owned dof should be in local range."));

  return 0;
}
