//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check CellAccessor::is_locally_owned

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


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_cube(tr, 3);

  typename Triangulation<dim,dim>::active_cell_iterator cell;

  for (cell = tr.begin_active();
       cell != tr.end();
       ++cell)
    {
      if (cell->is_locally_owned())
	{
	  if (myid==0)
	    deallog << cell << ": locally owned"
		    << std::endl;
	  Assert (!cell->is_ghost() && !cell->is_artificial(),
		  ExcInternalError());
	}
      else if (cell->is_ghost())
	{
	  if (myid==0)
	    deallog << cell << ": ghost"
		    << std::endl;
	  Assert (!cell->is_locally_owned() && !cell->is_artificial(),
		  ExcInternalError());
	}
      else if (cell->is_artificial())
	{
	  if (myid==0)
	    deallog << cell << ": artificial"
		    << std::endl;
	  Assert (!cell->is_locally_owned() && !cell->is_ghost(),
		  ExcInternalError());
	}
      else
	Assert (false, ExcInternalError());
    }
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("is_locally_owned").c_str());
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
