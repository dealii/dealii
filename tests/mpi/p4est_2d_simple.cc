//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// test p4est grid interaction with one cell in 2d

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/utilities.h>


#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD);

  if (true)
    {
      if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
	deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);
				       //tr.refine_global(1);

      Assert(tr.n_active_cells()==1, ExcInternalError());

      if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
	deallog << "subdomainid = " << tr.begin_active()->subdomain_id() << std::endl;

      if (myid == numproc-1)
	{
	  Assert( tr.begin_active()->subdomain_id()==(unsigned int)myid, ExcInternalError() );
	}
      else
	{
	  Assert( tr.begin_active()->subdomain_id()==types::artificial_subdomain_id, ExcInternalError() );
	}


      write_vtk (tr, "p4est_2d_simple", "1");
    }


  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_2d_simple").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();


#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
