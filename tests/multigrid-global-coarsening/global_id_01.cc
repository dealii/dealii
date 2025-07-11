// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test CellIDTranslator.

#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.templates.h>

#include <set>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const MPI_Comm comm)
{
  Triangulation<dim> basetria;
  GridGenerator::subdivided_hyper_cube(basetria, 4);
  basetria.refine_global(4);

  const auto determine_n_coarse_cells = [&comm](auto &tria) {
    types::coarse_cell_id n_coarse_cells = 0;

    for (auto cell : tria.active_cell_iterators())
      if (!cell->is_artificial())
        n_coarse_cells =
          std::max(n_coarse_cells, cell->id().get_coarse_cell_id());

    return Utilities::MPI::max(n_coarse_cells, comm) + 1;
  };

  // create translator: CellID <-> unique ID
  internal::CellIDTranslator<dim> cell_id_translator(
    determine_n_coarse_cells(basetria), basetria.n_global_levels());


  for (auto cell : basetria.cell_iterators())
    {
      Assert(cell->id() == cell_id_translator.to_cell_id(
                             cell_id_translator.translate(cell)),
             ExcNotImplemented());

      if (cell->has_children())
        {
          for (unsigned int c = 0; c < cell->n_children(); ++c)
            {
              Assert(cell->child(c)->id() ==
                       cell_id_translator.to_cell_id(
                         cell_id_translator.translate(cell, c)),
                     ExcNotImplemented());
            }
        }
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  test<2, 2>(comm);
}
