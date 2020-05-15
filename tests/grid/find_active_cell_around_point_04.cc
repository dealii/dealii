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

// Check for a bug discovered in step-70

#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize init(argc, argv);
  MPILogInitAll                    log;

  parallel::distributed::Triangulation<2> tria(
    MPI_COMM_WORLD,
    typename Triangulation<2>::MeshSmoothing(
      Triangulation<2>::smoothing_on_refinement |
      Triangulation<2>::smoothing_on_coarsening));

  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(5);
  GridTools::Cache<2> cache(tria);

  Point<2> p(0.239367, 0.341747);

  auto res1 = GridTools::find_active_cell_around_point(cache, p).first;
  auto res2 = GridTools::find_active_cell_around_point(tria, p);

  Assert(res1 == res2, ExcInternalError());
}
