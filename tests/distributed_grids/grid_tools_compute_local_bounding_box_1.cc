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
void test_hypercube(unsigned int ref, unsigned int max_boxes)
{
  const MPI_Comm &mpi_communicator = MPI_COMM_WORLD;
  deallog << "Testing hypercube for spacedim = " << spacedim << " ref " << ref << std::endl;

  parallel::distributed::Triangulation<spacedim> tria(mpi_communicator);
  GridGenerator::hyper_cube (tria);
  tria.refine_global(4);

  std::vector< BoundingBox<spacedim> > local_bbox = parallel::GridTools::compute_locally_owned_bounding_box<dim, spacedim>
                                                    ( tria, ref, max_boxes);
  if ( local_bbox.size() > 0)
    {
      deallog << "Computed Bounding Boxes:" << std::endl;
      for (BoundingBox<spacedim> b_box: local_bbox)
        {
          deallog << b_box.get_boundary_points().first << std::endl;
          deallog << b_box.get_boundary_points().second << std::endl;
          deallog << std::endl;
        }
    }
  else
    deallog << "Has no locally owned children cells" << std::endl;

  //Checking if all the points are inside the bounding boxes
  bool check = true;

  typename parallel::distributed::Triangulation< dim, spacedim >::active_cell_iterator
  cell = tria.begin_active();
  typename parallel::distributed::Triangulation< dim, spacedim >::active_cell_iterator
  endc = tria.last_active();

  //Looking if every point is at least inside a bounding box
  for (; cell<endc; ++cell)
    if (cell->is_locally_owned())
      for (unsigned int v=0; v<GeometryInfo<spacedim>::vertices_per_cell; ++v)
        {
          bool inside_a_box = false;
          for (unsigned int i=0; i< local_bbox.size(); ++i)
            {
              if (local_bbox[i].point_inside(cell->vertex(v)))
                inside_a_box = true;
            }
          if (! inside_a_box)
            {
              check = false;
              deallog << "Point outside " << cell->vertex(v) << std::endl;
              break;
            }
        }

  deallog << "Bounding Boxes surround locally_owned cells: " << check << std::endl;

  deallog << "Current test END"  << std::endl;
}

int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;

  unsigned int max_boxes = 4;

  test_hypercube<2> (0,max_boxes);
  test_hypercube<2> (1,max_boxes);
  test_hypercube<2> (2,max_boxes);
  test_hypercube<2> (3,max_boxes);
  test_hypercube<3> (1,max_boxes);
  test_hypercube<3> (2,max_boxes);
  test_hypercube<3> (3,max_boxes);


  return 0;
}
