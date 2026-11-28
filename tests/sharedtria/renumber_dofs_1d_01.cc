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

// Check DoFHandler::renumber_dofs() on a 1d parallel::shared::Triangulation
// with artificial cells. Such a triangulation uses the ParallelDistributed
// policy, which renumbers only the locally owned degrees of freedom. Reverse
// the locally owned degrees of freedom and check the indices on locally owned
// and ghost cells.

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


void
test()
{
  parallel::shared::Triangulation<1> tria(MPI_COMM_WORLD,
                                          Triangulation<1>::none,
                                          /*allow_artificial_cells=*/true);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(4);

  DoFHandler<1> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<1>(2));

  const IndexSet locally_owned = dof_handler.locally_owned_dofs();
  const types::global_dof_index n_locally_owned = locally_owned.n_elements();

  // The new index of a locally owned degree of freedom:
  const auto reversed = [&](const types::global_dof_index i) {
    return locally_owned.nth_index_in_set(n_locally_owned - 1 -
                                          locally_owned.index_within_set(i));
  };

  // Store the old indices on locally owned cells:
  std::map<std::pair<int, int>, std::vector<types::global_dof_index>>
    old_indices;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        auto &indices = old_indices[{cell->level(), cell->index()}];
        indices.resize(cell->get_fe().n_dofs_per_cell());
        cell->get_dof_indices(indices);
      }

  std::vector<types::global_dof_index> new_numbers(n_locally_owned);
  for (types::global_dof_index i = 0; i < n_locally_owned; ++i)
    new_numbers[i] = reversed(locally_owned.nth_index_in_set(i));
  dof_handler.renumber_dofs(new_numbers);

  AssertThrow(dof_handler.locally_owned_dofs() == locally_owned,
              ExcInternalError());

  // Locally owned degrees of freedom on locally owned cells have the
  // reversed index:
  for (const auto &[cell_id, indices] : old_indices)
    {
      const typename DoFHandler<1>::active_cell_iterator cell(&tria,
                                                              cell_id.first,
                                                              cell_id.second,
                                                              &dof_handler);
      std::vector<types::global_dof_index> new_indices(indices.size());
      cell->get_dof_indices(new_indices);
      for (unsigned int i = 0; i < indices.size(); ++i)
        if (locally_owned.is_element(indices[i]))
          AssertThrow(new_indices[i] == reversed(indices[i]),
                      ExcInternalError());
    }

  // All degrees of freedom on locally owned and ghost cells are valid and
  // locally relevant:
  const IndexSet locally_relevant =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  unsigned int n_ghost_cells = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (!cell->is_artificial())
      {
        n_ghost_cells += cell->is_ghost();
        std::vector<types::global_dof_index> indices(
          cell->get_fe().n_dofs_per_cell());
        cell->get_dof_indices(indices);
        for (const auto i : indices)
          AssertThrow(locally_relevant.is_element(i), ExcInternalError());
      }

  deallog << "n_dofs: " << dof_handler.n_dofs()
          << ", n_locally_owned_dofs: " << n_locally_owned
          << ", n_ghost_cells: " << n_ghost_cells << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test();
}
