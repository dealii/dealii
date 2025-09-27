// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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

  try
    {
      auto res1 = GridTools::find_active_cell_around_point(cache, p).first;
      auto res2 = GridTools::find_active_cell_around_point(tria, p);
      Assert(res1 == res2, ExcInternalError());
      deallog << "found on one processor" << std::endl;
    }
  catch (...)
    {}
}
