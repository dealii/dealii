// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Test repeated local mesh refinement (test contributed in PR #10092).

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include "../grid/tests.h"

using namespace dealii;


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  Triangulation<3, 3> tria;

  GridGenerator::subdivided_hyper_rectangle(tria,
                                            {{1.0, 1.0, 1.0, 1.0},
                                             {1.0, 1.0, 1.0, 1.0},
                                             {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}},
                                            Point<3>(-2.0, -2.0, -3.0),
                                            Point<3>(+2.0, +2.0, +3.0));

  for (auto cell : tria.active_cell_iterators())
    if (cell->center()[2] < 0.0)
      cell->set_material_id(1);
    else
      cell->set_material_id(0);

  for (unsigned int i = 0; i < 2; ++i)
    {
      for (auto cell : tria.active_cell_iterators())
        {
          if (cell->material_id() == 1)
            continue;

          for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
            if (cell->at_boundary(f) == false &&
                cell->neighbor(f)->material_id() == 1)
              cell->set_refine_flag();
        }

      tria.execute_coarsening_and_refinement();
    }

  parallel::fullydistributed::Triangulation<3> tria_pft(MPI_COMM_WORLD);

  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria, MPI_COMM_WORLD);

  tria_pft.create_triangulation(construction_data);

  deallog << "OK!" << std::endl;
}
