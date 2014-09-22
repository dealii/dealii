// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// refine bottom-left cell after one global refinement of a square in 2d and check p4est-output

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>


#include <fstream>

template <class TRIA>
void check (TRIA &tr)
{
  typename TRIA::cell_iterator cell = tr.begin(),
                               endc = tr.end();

  for (; cell!=endc; ++cell)
    {
      deallog << "cell level=" << cell->level() << " index=" << cell->index();
      if (!cell->has_children())
        deallog << " subdomain: " << cell->subdomain_id();
      deallog << std::endl;
    }

  deallog << "OK" << std::endl;
}


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (true)
    {
      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);
      tr.refine_global(1);
      tr.begin_active()->set_refine_flag();

      tr.execute_coarsening_and_refinement ();

      if (myid == 0)
        {
          deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
        }

      const unsigned int checksum = tr.get_checksum ();
      deallog << "Checksum: "
              << checksum
              << std::endl;

      check(tr);
    }



  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

}
