// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
#include "../tests.h"

#include "rotated_periodicity_common.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();
  for (const unsigned int angle : {0, 45, 90, 180})
    for (const unsigned int degree : {1, 2, 3})
      for (const unsigned int components : {1, 2, 3})
        {
          Triangulation<2> tria(
            Triangulation<2>::limit_level_difference_at_vertices);
          RotatedPeriodicityTest::check_mesh(
            tria, MPI_COMM_SELF, angle, degree, components);
        }
  Triangulation<3> tria(Triangulation<3>::limit_level_difference_at_vertices);
  RotatedPeriodicityTest::check_mesh(tria, MPI_COMM_SELF, 45, 2, 4);
  deallog << "Analytic periodic constraints and transfers: OK" << std::endl;
}
