// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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



// recursively refine a 2d mesh to a very high level

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <ostream>

#include "../tests.h"

template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (true)
    {
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);
      //      tr.refine_global(1);

      int level = 0;

      for (int i = 0; i < 29; ++i)
        {
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            deallog << "refine_loop... level=" << level << std::endl;

          if (myid == 0)
            tr.begin_active(level)->set_refine_flag();

          deallog.push("crap");
          tr.execute_coarsening_and_refinement();
          deallog.pop();


          if (myid == 0)
            {
              deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
            }

          level++;
        }
    }

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog.depth_file(3);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();
}
