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



// Test interaction with p4est with a few simple coarse grids in 2d

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
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/utilities.h>


#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (true)
    {
      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
      GridIn<dim> gi;
      gi.attach_triangulation (tr);
      std::ifstream in (SOURCE_DIR "/../deal.II/grid_in_02/2d.xda");
      try
        {
          gi.read_xda (in);
        }
      catch (const typename Triangulation<dim>::DistortedCellList &distorted_cells)
        {
          // ignore distorted cells
          deallog << distorted_cells.distorted_cells.size()
                  << " distorted cells after creating mesh."
                  << std::endl;
        }

      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        deallog << "subdomainid = "
                << tr.begin_active()->subdomain_id()
                << std::endl;

//      std::vector< unsigned int > cell_subd;
//      cell_subd.resize(tr.n_active_cells());

//      GridTools::get_subdomain_association(tr, cell_subd);
//       for (unsigned int i=0;i<tr.n_active_cells();++i)
//  deallog << cell_subd[i] << " ";
//       deallog << std::endl;

      if (myid == 0)
        {
          deallog << "#cells = " << tr.n_global_active_cells() << std::endl;

          Assert(tr.n_global_active_cells() == tr.n_active_cells(),
                 ExcInternalError() );
        }

      const unsigned int checksum = tr.get_checksum ();
      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        deallog << "Checksum: "
                << checksum
                << std::endl;
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

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();

}
