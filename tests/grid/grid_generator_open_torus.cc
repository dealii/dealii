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

// Short test to validate GridGenerator::torus<3,3>() with angle < 2*pi
// To test that the grid and the manifolds are generated correctly for
// open torus, we compare the aspect ratio of a full torus and an open
// torus for meshes with the same number of elements per arc length.

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"


void
test()
{
  Triangulation<3>   tria;
  const unsigned int n_cells_toroidal = 9;
  const double       angle            = 2.0 * numbers::PI;
  const double       R = 3., r = 1.;
  GridGenerator::torus(tria, R, r, n_cells_toroidal, angle);

  Triangulation<3>   tria_open;
  const unsigned int factor = 3;
  GridGenerator::torus(
    tria_open, R, r, n_cells_toroidal / factor, angle / (double)factor);

  const MappingQ<3> mapping(3);
  const QGauss<3>   gauss(4);

  const double ar_full_torus =
    GridTools::compute_maximum_aspect_ratio(mapping, tria, gauss);
  const double ar_open_torus =
    GridTools::compute_maximum_aspect_ratio(mapping, tria_open, gauss);

  deallog << "N_active_cells_full = " << tria.n_global_active_cells()
          << std::endl;
  deallog << "N_active_cells_open = " << tria_open.n_global_active_cells()
          << std::endl;
  deallog << "| AR_full - AR_open | = "
          << std::abs(ar_full_torus - ar_open_torus) << std::endl;
}

int
main(int argc, char **argv)
{
  initlog();
  deallog << std::scientific << std::setprecision(6);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

      test();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
