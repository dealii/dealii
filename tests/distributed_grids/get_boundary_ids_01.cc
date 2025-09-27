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


// Test parallel::TriangulationBase::get_boundary_ids()

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria, 4);

  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        {
          if (std::fabs(cell->face(face)->center()[0] - 0.0) < 1e-12)
            cell->face(face)->set_all_boundary_ids(1);
          if (std::fabs(cell->face(face)->center()[0] - 1.0) < 1e-12)
            cell->face(face)->set_all_boundary_ids(2);
          if (std::fabs(cell->face(face)->center()[1] - 0.0) < 1e-12)
            cell->face(face)->set_all_boundary_ids(3);
          if (std::fabs(cell->face(face)->center()[1] - 1.0) < 1e-12)
            cell->face(face)->set_all_boundary_ids(4);
        }

  for (const auto i : tria.get_boundary_ids())
    deallog << i << ' ';
  deallog << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2>();
}
