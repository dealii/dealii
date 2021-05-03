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
          if (std::fabs(cell->face(face)->center()(0) - 0.0) < 1e-12)
            cell->face(face)->set_all_boundary_ids(1);
          if (std::fabs(cell->face(face)->center()(0) - 1.0) < 1e-12)
            cell->face(face)->set_all_boundary_ids(2);
          if (std::fabs(cell->face(face)->center()(1) - 0.0) < 1e-12)
            cell->face(face)->set_all_boundary_ids(3);
          if (std::fabs(cell->face(face)->center()(1) - 1.0) < 1e-12)
            cell->face(face)->set_all_boundary_ids(4);
        }

  for (const auto i : tria.get_boundary_ids())
    deallog << i << " ";
  deallog << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2>();
}
