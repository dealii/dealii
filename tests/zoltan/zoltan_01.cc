// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check GridTools::partition_triangulation using ZOLTAN as partitioner
// Test 1 (metis_01.cc) of metis is used as model for this test

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>


template <int dim>
void
test ()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (5-dim);

  // subdivides into 5 subdomains
  deallog << "Partitioning" << std::endl;
  GridTools::partition_triangulation (5, triangulation, SparsityTools::Partitioner::zoltan);
  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    deallog << cell << ' ' << cell->subdomain_id() << std::endl;

}



int
main (int argc, char **argv)
{
  //Initialize MPI and Zoltan
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);
  MPILogInitAll all;

  //tests
  test<1> ();
  test<2> ();
  test<3> ();
}
