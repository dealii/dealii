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



// refine bottom-left cell after one global refinement of a square in 2d and
// check p4est-output

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <class TRIA>
void
check(TRIA &tr)
{
  typename TRIA::cell_iterator cell = tr.begin(), endc = tr.end();

  for (; cell != endc; ++cell)
    {
      deallog << "cell level=" << cell->level() << " index=" << cell->index();
      if (!cell->has_children())
        deallog << " subdomain: " << cell->subdomain_id();
      deallog << std::endl;
    }

  deallog << "OK" << std::endl;
}


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
      tr.refine_global(1);
      tr.begin_active()->set_refine_flag();

      tr.execute_coarsening_and_refinement();

      if (myid == 0)
        {
          deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
        }

      const unsigned int checksum = tr.get_checksum();
      deallog << "Checksum: " << checksum << std::endl;

      check(tr);
    }



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
