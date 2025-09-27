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



// check CellAccessor::is_locally_owned

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid     = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_cube(tr, 3);

  typename Triangulation<dim, dim>::active_cell_iterator cell;

  for (cell = tr.begin_active(); cell != tr.end(); ++cell)
    {
      if (cell->is_locally_owned())
        {
          if (myid == 0)
            deallog << cell << ": locally owned" << std::endl;
          Assert(!cell->is_ghost() && !cell->is_artificial(),
                 ExcInternalError());
        }
      else if (cell->is_ghost())
        {
          if (myid == 0)
            deallog << cell << ": ghost" << std::endl;
          Assert(!cell->is_locally_owned() && !cell->is_artificial(),
                 ExcInternalError());
        }
      else if (cell->is_artificial())
        {
          if (myid == 0)
            deallog << cell << ": artificial" << std::endl;
          Assert(!cell->is_locally_owned() && !cell->is_ghost(),
                 ExcInternalError());
        }
      else
        Assert(false, ExcInternalError());
    }
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

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();
}
