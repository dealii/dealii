// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// A testcase by Andrea Bonito (in modified form): Given a particular
// set of refinement indicators, we get no refinement at all because
// the threshold is computed to be *larger* than the largest
// refinement indicator

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // create a mesh with 120 cells (because that's what Andrea's
  // original testcase had)
  parallel::distributed::Triangulation<2> triangulation(MPI_COMM_WORLD);
  std::vector<unsigned int>               subdivisions(2);
  subdivisions[0] = 120;
  subdivisions[1] = 1;
  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            subdivisions,
                                            Point<2>(),
                                            Point<2>(1, 1));
  // initialize the refinement indicators with a set of particular values from
  // the original testcase
  const double values[] = {
    1.48589e-08, 3.31859e-06, 3.31859e-06, 0.0025918,   3.31859e-06,
    0.0025918,   0.0025918,   3.31859e-06, 1.48589e-08, 0.0025918,
    3.31859e-06, 0.0025918,   3.31859e-06, 0.0025918,   3.31859e-06,
    0.0025918,   1.48589e-08, 3.31859e-06, 0.0025918,   3.31859e-06,
    0.0025918,   0.0025918,   3.31859e-06, 3.31859e-06, 1.48589e-08,
    0.0025918,   0.0025918,   3.31859e-06, 3.31859e-06, 0.0025918,
    0.0025918,   1.48589e-08, 3.31859e-06, 3.31859e-06, 0.0025918,
    0.0025918,   3.31859e-06, 0.0025918,   3.31859e-06, 1.48589e-08,
    0.0025918,   3.31859e-06, 0.0025918,   3.31859e-06, 0.0025918,
    3.31859e-06, 0.0025918,   1.48589e-08, 3.31859e-06, 0.0025918,
    0.0025918,   3.31859e-06, 0.0025918,   3.31859e-06, 3.31859e-06,
    1.48589e-08, 0.0120662,   0.0446999,   0.0446999,   0.0644361,
    0.0446999,   0.0644361,   0.0644361,   3.14051,     0.0446999,
    0.0120662,   0.0644361,   0.0446999,   0.0644361,   0.0446999,
    3.14051,     0.0644361,   0.0446999,   0.0644361,   0.0120662,
    0.0446999,   0.0644361,   3.14051,     0.0446999,   0.0644361,
    0.0644361,   0.0446999,   0.0446999,   0.0120662,   3.14051,
    0.0644361,   0.0644361,   0.0446999,   0.0446999,   0.0644361,
    0.0644361,   3.14051,     0.0120662,   0.0446999,   0.0446999,
    0.0644361,   0.0644361,   0.0446999,   3.14051,     0.0644361,
    0.0446999,   0.0120662,   0.0644361,   0.0446999,   0.0644361,
    3.14051,     0.0446999,   0.0644361,   0.0446999,   0.0644361,
    0.0120662,   0.0446999,   3.14051,     0.0644361,   0.0644361,
    0.0446999,   0.0644361,   0.0446999,   0.0446999,   0.0120662};

  Vector<double> estimated_error_per_cell(&values[0], &values[0] + 120);

  deallog << "Number of cells before refinement: "
          << triangulation.n_active_cells() << std::endl;

  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
    triangulation, estimated_error_per_cell, 0.3, 0.0);
  triangulation.execute_coarsening_and_refinement();

  deallog << "Number of cells after refinement: "
          << triangulation.n_active_cells() << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      test();
    }
  else
    test();
}
