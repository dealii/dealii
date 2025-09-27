// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check GridTools::partition_triangulation using ZOLTAN as partitioner
// Test 1 (metis_01.cc) of metis is used as model for this test

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(5 - dim);

  // subdivides into 5 subdomains
  deallog << "Partitioning" << std::endl;
  GridTools::partition_triangulation(5,
                                     triangulation,
                                     SparsityTools::Partitioner::zoltan);
  for (typename Triangulation<dim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    deallog << cell << ' ' << cell->subdomain_id() << std::endl;
}



int
main(int argc, char **argv)
{
  // Initialize MPI and Zoltan
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  MPILogInitAll                    all;

  // tests
  test<1>();
  test<2>();
  test<3>();
}
