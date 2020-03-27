// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"


void
test()
{
  Triangulation<3>   tria;
  unsigned int const n_cells_toroidal = 9;
  double const       angle            = 2.0 * numbers::PI;
  double const       R = 3., r = 1.;
  GridGenerator::torus(tria, R, r, n_cells_toroidal, angle);

  Triangulation<3>   tria_open;
  unsigned int const factor = 3;
  GridGenerator::torus(
    tria_open, R, r, n_cells_toroidal / factor, angle / (double)factor);

  MappingQGeneric<3> const mapping(3);
  QGauss<3> const          gauss(4);

  double const ar_full_torus =
    GridTools::compute_maximum_aspect_ratio(tria, mapping, gauss);
  double const ar_open_torus =
    GridTools::compute_maximum_aspect_ratio(tria_open, mapping, gauss);

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
  catch (std::exception &exc)
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
