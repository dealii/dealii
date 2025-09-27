// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test that Cache::get_covering_rtree() works

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/boost_adaptors/bounding_box.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"

namespace bgi = boost::geometry::index;

template <int dim, int spacedim = dim>
void
test(unsigned int ref)
{
  const MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  parallel::distributed::Triangulation<dim, spacedim> tria(mpi_communicator);
  GridGenerator::hyper_ball(tria);
  tria.refine_global(ref);

  GridTools::Cache<dim, spacedim> cache(tria);

  const auto &tree  = cache.get_covering_rtree();
  const auto  point = random_point<dim>();

  deallog << "Query point: " << point << std::endl;

  for (const auto p : tree | bgi::adaptors::queried(bgi::intersects(point)))
    deallog << "Processor " << p.second << " may own the point " << point
            << ": it is within " << p.first.get_boundary_points().first << "; "
            << p.first.get_boundary_points().second << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>(3);

  return 0;
}
