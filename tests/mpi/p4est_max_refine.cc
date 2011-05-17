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


// recursively refine a 2d mesh to a very high level

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
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);

  if (true)
    {
      if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
	deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);
//      tr.refine_global(1);

      int level = 0;
      
      for (int i=0;i<29;++i)
	{
	  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
	    deallog << "refine_loop... level=" << level  << std::endl;
	  
	  if (myid==0)
	    tr.begin_active(level)->set_refine_flag();

	  deallog.push("crap");
	  tr.execute_coarsening_and_refinement ();
	  deallog.pop();
	  

	  if (myid == 0)
	    {
	      deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
	    }

	  level++;
//	  write_vtk (tr, "p4est_max_refine", Utilities::int_to_string(i,4).c_str());

	}
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


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_max_refine").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.depth_file(3);
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
