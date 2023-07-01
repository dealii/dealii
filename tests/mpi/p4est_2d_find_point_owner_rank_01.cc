// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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
   * hyper_L. We want to find the mpi rank of a single fixed point. On all
   * processes we must find the same owner.
   */
  const std::vector<Point<2>> point(1, Point<2>(0.23, 0.77));

  deallog.push("2D-Test");
  test<2>(point);
  deallog.pop();
}
