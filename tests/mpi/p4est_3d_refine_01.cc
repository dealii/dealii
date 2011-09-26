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
	      const unsigned int x = rand() % flags.size();
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
      std::ofstream logfile(output_file_for_mpi("p4est_3d_refine_01").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    test<3>();


#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
