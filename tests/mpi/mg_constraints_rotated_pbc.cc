// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
#include "../tests.h"

#include "../multigrid/rotated_periodicity_common.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  mpi_initlog();
  for (const unsigned int angle : {0, 45, 90, 180})
    for (const unsigned int degree : {1, 2, 3})
      for (const unsigned int components : {1, 2, 3})
        {
          parallel::distributed::Triangulation<2> tria(
            MPI_COMM_WORLD,
            Triangulation<2>::limit_level_difference_at_vertices,
            parallel::distributed::Triangulation<
              2>::construct_multigrid_hierarchy);
          RotatedPeriodicityTest::check_mesh(
            tria, MPI_COMM_WORLD, angle, degree, components);
        }
  parallel::distributed::Triangulation<3> tria(
    MPI_COMM_WORLD,
    Triangulation<3>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<3>::construct_multigrid_hierarchy);
  RotatedPeriodicityTest::check_mesh(tria, MPI_COMM_WORLD, 45, 2, 4);
  deallog << "Analytic periodic constraints and transfers: OK" << std::endl;
}
