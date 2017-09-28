// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test for the function compute_locally_owned_bounding_box : on various domains and
// mpi configurations to vary the shape of the various domains


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/distributed/grid_tools.h>
#include <deal.II/distributed/grid_refinement.h>

template <int dim, int spacedim=dim >
void test_hypercube(unsigned int ref)
{
  const MPI_Comm &mpi_communicator = MPI_COMM_WORLD;
  deallog << "Testing hypercube for spacedim = " << spacedim << " ref " << ref << std::endl;

  parallel::distributed::Triangulation<spacedim> tria(mpi_communicator);
  GridGenerator::hyper_cube (tria);
  tria.refine_global(4);

  for (typename dealii::parallel::Triangulation<dim, spacedim>::cell_iterator cell: tria.cell_iterators_on_level(ref))
    {
      bool has_locally_owned = true;
      BoundingBox<spacedim> local_bbox = parallel::GridTools::compute_cell_locally_owned_bounding_box<dim, spacedim>(cell,has_locally_owned);

      if (has_locally_owned)
        {
          deallog << "Computed Bounding Box:" << std::endl;
          deallog << local_bbox.get_boundary_points().first << std::endl;
          deallog << local_bbox.get_boundary_points().second << std::endl;
        }
      else
        deallog << "Has no locally owned children cells" << std::endl;
    }

  deallog << "END test hypercube for spacedim = " << spacedim << " ref " << ref << std::endl;
}

int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;

  test_hypercube<2> (0);
  test_hypercube<2> (1);
  test_hypercube<2> (2);
  test_hypercube<3> (0);
  test_hypercube<3> (1);
  test_hypercube<3> (2);

  return 0;
}
