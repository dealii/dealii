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

// Check that DoFAccessor::vertex_dof_index() returns
// numbers::invalid_dof_index on vertices that are only adjacent to artificial
// cells, and a valid index on all other vertices. Check this with and without
// hp capabilities.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


template <int dim>
void
test(const hp::FECollection<dim> &fe)
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  std::vector<unsigned int> repetitions(dim, 1);
  repetitions[0] = 8;
  Point<dim> top_right;
  for (unsigned int d = 0; d < dim; ++d)
    top_right[d] = (d == 0 ? 8. : 1.);
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            repetitions,
                                            Point<dim>(),
                                            top_right);
  tria.refine_global(1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  std::vector<bool> vertex_has_storage(tria.n_vertices(), false);
  for (const auto &cell : tria.active_cell_iterators())
    if (!cell->is_artificial())
      for (const auto v : cell->vertex_indices())
        vertex_has_storage[cell->vertex_index(v)] = true;

  unsigned int n_invalid = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_artificial())
      for (const auto v : cell->vertex_indices())
        {
          const bool is_valid =
            cell->vertex_dof_index(v, 0, 0) != numbers::invalid_dof_index;
          AssertThrow(is_valid == vertex_has_storage[cell->vertex_index(v)],
                      ExcInternalError());
          n_invalid += !is_valid;
        }

  deallog << "hp: " << dof_handler.has_hp_capabilities()
          << ", invalid vertex DoF indices on artificial cells: " << n_invalid
          << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>(hp::FECollection<2>(FE_Q<2>(1)));
  test<2>(hp::FECollection<2>(FE_Q<2>(1), FE_Q<2>(2)));
  deallog.pop();

  deallog.push("3d");
  test<3>(hp::FECollection<3>(FE_Q<3>(1)));
  test<3>(hp::FECollection<3>(FE_Q<3>(1), FE_Q<3>(2)));
  deallog.pop();
}
