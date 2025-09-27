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



// test n_levels and iterators on those levels
// begin(i) should work even if we have no cells on that level

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  if (myid == 0)
    deallog << "#cells = " << tr.n_global_active_cells() << std::endl;

  deallog << "proc " << myid << ", n levels " << tr.n_levels() << std::endl;
  deallog << "proc " << myid << ", n global levels " << tr.n_global_levels()
          << std::endl;

  deallog << "begin().cell_index " << tr.begin()->index() << std::endl;
  deallog << "begin(0).cell_index " << tr.begin(0)->index() << std::endl;

  deallog << "begin(1)==end(1)? " << (tr.begin(1) == tr.end(1)) << std::endl;


  deallog << "subdomainid = " << tr.begin_active()->subdomain_id() << std::endl;

  // if (myid!=0)
  //   Assert(tr.begin(1)==tr.end(1), ExcInternalError());

  const unsigned int checksum = tr.get_checksum();
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "Checksum: " << checksum << std::endl;

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();
}
