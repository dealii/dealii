// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Test CellIDTranslator.

#include <deal.II/grid/cell_id_translator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const MPI_Comm comm)
{
  Triangulation<dim> basetria;
  GridGenerator::subdivided_hyper_cube(basetria, 4);
  basetria.refine_global(4);

  // create translator: CellID <-> unique ID
  internal::CellIDTranslator<dim> cell_id_translator(basetria);
  for (const auto &cell : basetria.cell_iterators())
    {
      Assert(cell->id() == cell_id_translator.to_cell_id(
                             cell_id_translator.translate(cell)),
             ExcInternalError());

      if (cell->has_children())
        {
          for (unsigned int c = 0; c < cell->n_children(); ++c)
            {
              Assert(cell->child(c)->id() ==
                       cell_id_translator.to_cell_id(
                         cell_id_translator.translate(cell, c)),
                     ExcInternalError());
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
