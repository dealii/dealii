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



// recursively refine a 3d mesh

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
#include <ostream>

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

      while (tr.n_active_cells() < 20000/Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
        {
          if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
            deallog << "refine_loop..." << std::endl;
          std::vector<bool> flags (tr.n_active_cells(), false);

          // refine one fifth of all cells each
          // time (but at least one).
          // note that only the own marked cells
          // will be refined.
          for (unsigned int i=0; i<tr.n_active_cells() / 5 + 1; ++i)
            {
              const unsigned int x = Testing::rand() % flags.size();
              flags[x] = true;
            }

          unsigned int index=0;
          unsigned int locals=0;

          for (typename Triangulation<dim>::active_cell_iterator
               cell = tr.begin_active();
               cell != tr.end(); ++cell, ++index)
            if (flags[index])
              {
                if (cell->subdomain_id()==myid)
                  ++locals;
                cell->set_refine_flag();
              }

          //and atleast one cell:
          if (!locals)
            {

              for (typename Triangulation<dim>::active_cell_iterator
                   cell = tr.begin_active();
                   cell != tr.end(); ++cell)
                if (cell->subdomain_id()==myid)
                  {
                    cell->set_refine_flag();
                    ++locals;
                    if (locals>5)


                      break;
                  }
            }




          Assert (index == tr.n_active_cells(), ExcInternalError());
          tr.execute_coarsening_and_refinement ();

          unsigned int checksum = tr.get_checksum ();
          if (myid == 0)
            {
              deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
              deallog << "Checksum: "
                      << checksum
                      << std::endl;
            }
        }
    }

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    test<3>();

}
