// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check GridTools::Cache::get_locally_owned_cell_bounding_boxes_rtree()


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
test()
{
  deallog << "Testing for dim = " << dim << std::endl;

  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::Settings::partition_zoltan);
  GridGenerator::hyper_cube(tria, -3, -2);
  tria.refine_global(std::max(6 - dim, 3));


  GridTools::Cache<dim> cache(tria);

  const auto &global_description =
    cache.get_locally_owned_cell_bounding_boxes_rtree();

  // Extract from the cache a list of all bounding boxes with process owners:
  std::vector<std::pair<BoundingBox<dim>,
                        typename Triangulation<dim>::active_cell_iterator>>
    test_results;
  global_description.query(bgi::satisfies([](const auto &) { return true; }),
                           std::back_inserter(test_results));

  // Loop over all bounding boxes and output them. We expect this list
  // to be identical on all processes, but do not actually check this
  // (other than by the fact that the output file records this).
  for (const auto &bb_and_owner : test_results)
    {
      const auto &bd_points = bb_and_owner.first.get_boundary_points();
      deallog << "  Bounding box: p1 " << bd_points.first << " p2 "
              << bd_points.second << " owning cell: " << bb_and_owner.second
              << std::endl;
    }
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

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
