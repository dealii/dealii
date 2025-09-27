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



// test p4est grid interaction with one cell in 2d

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (true)
    {
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);
      // tr.refine_global(1);

      Assert(tr.n_active_cells() == 1, ExcInternalError());

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "subdomainid = " << tr.begin_active()->subdomain_id()
                << std::endl;

      if (myid == numproc - 1)
        {
          Assert(tr.begin_active()->subdomain_id() == (unsigned int)myid,
                 ExcInternalError());
        }
      else
        {
          Assert(tr.begin_active()->subdomain_id() ==
                   numbers::artificial_subdomain_id,
                 ExcInternalError());
        }


      const unsigned int checksum = tr.get_checksum();
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "Checksum: " << checksum << std::endl;
    }


  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

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
