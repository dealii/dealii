// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * Test MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence.
 */

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, -1.0, +1.0);

  tria.refine_global(1);

  // refine cells in the first quadrant
  for (unsigned int i = 1; i < 3; ++i)
    {
      for (auto cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            bool flag = true;
            for (int d = 0; d < dim; ++d)
              if (cell->center()[d] > 0.0)
                flag = false;
            if (flag)
              cell->set_refine_flag();
          }
      tria.execute_coarsening_and_refinement();
    }

  const auto trias =
    MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(tria);

  for (const auto &tria : trias)
    {
      for (const auto &cell : tria->active_cell_iterators())
        if (cell->is_locally_owned())
          deallog << cell->id() << std::endl;
      deallog << std::endl;
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  test<2>();
}
