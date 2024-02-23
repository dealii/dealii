// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <vector>

#include "../tests.h"

template <int dim>
void
test(MPI_Comm comm)
{
  deallog << "dim = " << dim << std::endl;
  parallel::shared::Triangulation<dim> tria(comm);

  std::vector<typename GridTools::PeriodicFacePair<
    typename parallel::shared::Triangulation<dim>::cell_iterator>>
    matched_pairs;

  unsigned int              num_refinements = 1 << 4;
  Point<dim>                p1;
  Point<dim>                p2;
  std::vector<unsigned int> repetitions(dim);

  for (unsigned int i = 0; i < dim; ++i)
    {
      p1[i]          = 0.;
      p2[i]          = 1.;
      repetitions[i] = num_refinements;
    }
  GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2, true);
  // Collect periodic faces in the x-direction
  GridTools::collect_periodic_faces(tria, 0, 1, 0, matched_pairs);
  // Check the size of the matched_pairs vector
  deallog << matched_pairs.size() << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();
  MPI_Comm comm = MPI_COMM_WORLD;
  test<1>(comm);
  test<2>(comm);
  test<3>(comm);

  return 0;
}
