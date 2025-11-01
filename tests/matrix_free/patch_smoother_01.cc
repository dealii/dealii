// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// check basic functionality of PatchStorage: initialization on a mesh with a
// single patch

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_patch_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/patch_storage.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // 1. Generate triangulation
  Triangulation<dim> triangulation(
    Triangulation<dim>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1); // Refine once

  // 2. Initialize FE, DoFHandler, Mapping
  const unsigned int fe_degree = 1;
  FE_DGQ<dim>        fe(fe_degree);
  DoFHandler<dim>    dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  const unsigned int n_lanes = VectorizedArray<double>::size();

  MappingQ1<dim> mapping; // Using Q1 mapping as in step-94

  // 3. Initialize constraints (needed for MatrixFree)
  AffineConstraints<double> constraints;
  IndexSet                  locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
  constraints.reinit(locally_relevant_dofs);
  // No boundary conditions or hanging nodes needed for this simple test yet
  constraints.close();

  // 4. Initialize MatrixFree for level 1
  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mg_level          = 1; // Target level 1
  additional_data.store_ghost_cells = true;

  std::shared_ptr<MatrixFree<dim, double>> mf_level_storage =
    std::make_shared<MatrixFree<dim, double>>();
  mf_level_storage->reinit(mapping,
                           dof_handler,
                           constraints,
                           QGauss<1>(fe_degree + 1),
                           additional_data);

  // 5. Build PatchStorage on level 1
  PatchStorage<MatrixFree<dim, double>> patch_storage(mf_level_storage);
  patch_storage.initialize();


  // Checks: there should only one patch

  Assert(patch_storage.n_patches() > 0, ExcInternalError());
  deallog << "  Level 1 PatchStorage initialized with "
          << patch_storage.n_patches() << " patches." << std::endl;


  // with cell ordered lexicographically, that is the order of childred
  auto       patch_cells = patch_storage.get_regular_patch(0).get_cells();
  const auto parent_cell = triangulation.begin(0);
  // Get the map from cell iterators (level 1) to MatrixFree indices

  deallog << "  Checking children of first coarse cell (" << parent_cell->id()
          << ") against first patch cells:" << std::endl;
  // Iterate through the children of the parent cell
  for (unsigned int child_idx = 0; child_idx < parent_cell->n_children();
       ++child_idx)
    {
      // Get the iterator for the child cell (which should be on level 1)
      const auto child_cell = parent_cell->child(child_idx);
      Assert(child_cell->level() == 1,
             ExcInternalError("Child cell is not on level 1"));
      Assert(child_cell->is_active(),
             ExcInternalError(
               "Child cell is not active")); // Cells on target level should be
                                             // active

      auto mf_cell_iterator =
        mf_level_storage->get_cell_iterator(patch_cells[child_idx] / n_lanes,
                                            patch_cells[child_idx] % n_lanes);
      Assert(
        mf_cell_iterator->index() == child_cell->index(),
        ExcMessage(
          "Cell index mismatch -- most likely patch construction went wrong."));
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  initlog();

  deallog << "Running test in 2D..." << std::endl;
  test<2>();

  deallog << "Running test in 3D..." << std::endl;
  test<3>();

  deallog << "Tests finished." << std::endl;

  return 0; // Indicate success
}
