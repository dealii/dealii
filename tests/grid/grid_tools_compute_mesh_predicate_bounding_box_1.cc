// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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


// test for the function compute_mesh_predicate_bounding_box : as a predicate
// is_locally_owned is used and, to vary the shapes, various mpi configurations
// are used on a distributed::parallel::triangulations


#include <deal.II/base/bounding_box.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"

// predicate: test if the cell is locally owned
template <int dim, int spacedim>
bool
pred_locally_owned(const typename parallel::distributed::
                     Triangulation<dim, spacedim>::cell_iterator &cell)
{
  return cell->is_locally_owned();
}

template <int dim, int spacedim = dim>
void
test_hypercube(unsigned int ref, unsigned int max_bbox)
{
  const MPI_Comm &mpi_communicator = MPI_COMM_WORLD;
  deallog << "Testing hypercube for spacedim = " << spacedim
          << " refinement: " << ref << " max number of boxes: " << max_bbox
          << std::endl;

  parallel::distributed::Triangulation<dim, spacedim> tria(mpi_communicator);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(4);

  typedef typename parallel::distributed::Triangulation<dim, spacedim>::
    active_cell_iterator cell_iterator;

  std::function<bool(const cell_iterator &)> predicate =
    pred_locally_owned<dim, spacedim>;

  std::vector<BoundingBox<spacedim>> local_bbox =
    GridTools::compute_mesh_predicate_bounding_box<
      parallel::distributed::Triangulation<dim, spacedim>>(
      tria, predicate, ref, true, max_bbox);

  if (local_bbox.size() > 0)
    {
      deallog << "Computed Bounding Boxes:" << std::endl;
      for (BoundingBox<spacedim> b_box : local_bbox)
        {
          deallog << b_box.get_boundary_points().first << std::endl;
          deallog << b_box.get_boundary_points().second << std::endl;
          deallog << std::endl;
        }
    }
  else
    deallog << "Has no locally owned children cells" << std::endl;

  // Checking if all the points are inside the bounding boxes
  bool check = true;

  typename parallel::distributed::Triangulation<dim,
                                                spacedim>::active_cell_iterator
    cell = tria.begin_active();
  typename parallel::distributed::Triangulation<dim,
                                                spacedim>::active_cell_iterator
    endc = tria.last_active();

  // Looking if every point is at least inside a bounding box
  for (; cell < endc; ++cell)
    if (cell->is_locally_owned())
      for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
        {
          bool inside_a_box = false;
          for (unsigned int i = 0; i < local_bbox.size(); ++i)
            {
              if (local_bbox[i].point_inside(cell->vertex(v)))
                inside_a_box = true;
            }
          if (!inside_a_box)
            {
              check = false;
              deallog << "Point outside " << cell->vertex(v) << std::endl;
              break;
            }
        }
  deallog << "Bounding Boxes surround locally_owned cells: " << check
          << std::endl;

  deallog << "Current test END" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test_hypercube<2>(0, 4);
  test_hypercube<2>(1, 4);
  test_hypercube<2>(2, 4);
  test_hypercube<2>(2, 2);
  // test_hypercube<2, 3>(0, 4);
  // test_hypercube<2, 3>(1, 4);
  // test_hypercube<2, 3>(2, 4);
  // test_hypercube<2, 3>(2, 2);
  test_hypercube<3>(0, 4);
  test_hypercube<3>(1, 4);
  test_hypercube<3>(2, 4);
  test_hypercube<3>(2, 2);

  return 0;
}
