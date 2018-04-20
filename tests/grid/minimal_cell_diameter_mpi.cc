// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>


template <int dim>
void test1 ()
{
  // test 1: hypercube
  if (true)
    {
      parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
      GridGenerator::hyper_cube(tria);

      for (unsigned int i=0; i<2; ++i)
        {
          tria.refine_global(2);
          deallog << dim << "d, "
                  << "min diameter: "
                  << GridTools::minimal_cell_diameter (tria)
                  << std::endl;
        };
    };

  // test 2: hyperball
  if (dim >= 2)
    {
      parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
      GridGenerator::hyper_ball(tria, Point<dim>(), 1);
      tria.reset_manifold(0);

      for (unsigned int i=0; i<2; ++i)
        {
          tria.refine_global(2);
          deallog << dim << "d, "
                  << "min diameter: "
                  << GridTools::minimal_cell_diameter (tria)
                  << std::endl;
        };
    };
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  MPILogInitAll mpi_init_log;

  test1<2> ();
  test1<3> ();

  return 0;
}

