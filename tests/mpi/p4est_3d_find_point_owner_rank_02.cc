// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/mpi.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include "../tests.h"

// contains the actual templated tests
#include "p4est_find_point_owner_rank.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;


  /*
   * We create a distributed mesh with three cells (tree roots), e.g., a
   * hyper_L. We want to find the mpi rank of a set of fixed points. On all
   * processes we must find the same owner ranks.
   */
  std::vector<Point<3>> points;
  points.emplace_back(Point<3>(0.23, 0.23, 0.23));
  points.emplace_back(Point<3>(0.23, 0.77, 0.23));
  points.emplace_back(Point<3>(0.77, 0.23, 0.23));
  points.emplace_back(Point<3>(0.77, 0.77, 0.23));
  points.emplace_back(Point<3>(0.23, 0.23, 0.77));
  points.emplace_back(Point<3>(0.23, 0.77, 0.77));
  points.emplace_back(Point<3>(0.77, 0.23, 0.77));
  points.emplace_back(Point<3>(0.77, 0.77, 0.77));

  // this point is outside the mesh
  points.emplace_back(Point<3>(1.1, 1.1, 1.1));

  deallog.push("3D-Test");
  test<3>(points);
  deallog.pop();
}
