// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
