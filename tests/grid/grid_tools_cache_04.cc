// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2015 by the deal.II authors
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

// Check extract used_vertices and find_closest_vertex using a cache

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/fe/mapping_q1.h>
#include <deal.II/base/mpi.templates.h>

// Testing GridTools::Cache with bounding boxes

template <int dim, int spacedim=dim >
void test_cache(unsigned int ref, bool allow_merge, unsigned int max_bbox)
{
  const MPI_Comm &comm = MPI_COMM_WORLD;
  deallog << "Testing GridTools::Cache for spacedim = " << spacedim << " and dimension " << dim << std::endl;

  parallel::distributed::Triangulation<spacedim> tria(comm);
  GridGenerator::hyper_cube (tria);
  tria.refine_global(4);

  GridTools::Cache<dim,spacedim> cache(tria,StaticMappingQ1<dim,spacedim>::mapping,comm,true,ref,allow_merge,max_bbox);

  auto global_bboxes = cache.get_global_bounding_boxes();

  const auto my_proc = dealii::Utilities::MPI::this_mpi_process(comm);

  if ( global_bboxes.size() > 0)
    {
      deallog << "Computed Bounding Boxes:" << std::endl;
      for (auto b_box_rk: global_bboxes)
        for (auto b_box: b_box_rk)
          {
            deallog << b_box.get_boundary_points().first << std::endl;
            deallog << b_box.get_boundary_points().second << std::endl;
            deallog << std::endl;
          }
    }
  else
    deallog << "Error: no global bounding boxes" << std::endl;

  deallog << "Test finished" << std::endl;
};


int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll init;

  test_cache<2,2> (1,true,4);

  return 0;
}

