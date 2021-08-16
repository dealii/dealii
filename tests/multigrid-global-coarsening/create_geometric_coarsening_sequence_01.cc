// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


/**
 * Test MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence.
 */

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"

using namespace dealii;

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
