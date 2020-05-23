// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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

// Validate grid_tools_cache for the creation of build global rtree


#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/signals2.hpp>

#include "../tests.h"


namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;


template <int dim>
void
test(unsigned int n_points)
{
  deallog << "Testing for dim = " << dim << std::endl;

  // Creating a grid in the square [0,1]x[0,1]
  parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(std::max(8 - dim, 3));


  GridTools::Cache<dim> cache(tria);

  const auto &global_description = cache.get_covering_rtree();

  // Creating the random points
  std::vector<Point<dim>> points;

  for (size_t i = 0; i < n_points; ++i)
    points.push_back(random_point<dim>());



  const auto cell_qpoint_map =
    GridTools::compute_point_locations(cache, points);
  const auto & cells   = std::get<0>(cell_qpoint_map);
  const auto & maps    = std::get<2>(cell_qpoint_map);
  unsigned int n_cells = cells.size();

  for (unsigned int c = 0; c < n_cells; ++c)
    {
      // We know the owner as we're working on a shared triangulation
      unsigned int cell_owner = cells[c]->subdomain_id();
      for (unsigned int idx : maps[c])
        {
          std::vector<std::pair<BoundingBox<dim>, unsigned int>> test_results;
          global_description.query(bgi::intersects(points[idx]),
                                   std::back_inserter(test_results));
          bool found = false;

          // Printing and checking function output
          for (const auto &ans : test_results)
            {
              const auto &bd_points = ans.first.get_boundary_points();
              deallog << "Bounding box: p1 " << bd_points.first << " p2 "
                      << bd_points.second << " rank owner: " << ans.second
                      << std::endl;
              // Check if the guess made from the tree contains the answer
              if (ans.second == cell_owner)
                found = true;
            }

          if (!found)
            deallog << "Error point " << points[idx] << " was not found."
                    << std::endl;
        }
    }
  deallog << "Test for dim " << dim << " finished." << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
#ifdef DEAL_II_WITH_MPI
  MPILogInitAll log;
#else
  initlog();
  deallog.push("0");
#endif

  test<1>(50);
  test<2>(80);
  test<3>(101);

  return 0;
}
