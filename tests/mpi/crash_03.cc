//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// extracted from count_dofs_per_block that crashed at one point


#include "../tests.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <distributed/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/intergrid_map.h>
#include <base/utilities.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>

#include <fstream>
#include <numeric>
#include <cstdlib>


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim>
    triangulation (MPI_COMM_WORLD,
		   Triangulation<dim>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global (2);

  for (unsigned int i=0; i<2; ++i)
    {
				       // refine one-fifth of cells randomly
      std::vector<bool> flags (triangulation.n_active_cells(), false);
      for (unsigned int k=0; k<flags.size()/5 + 1; ++k)
	flags[rand() % flags.size()] = true;
				       // make sure there's at least one that
				       // will be refined
      flags[0] = true;

				       // refine triangulation
      unsigned int index=0;
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = triangulation.begin_active();
	   cell != triangulation.end(); ++cell, ++index)
	if (flags[index])
	  cell->set_refine_flag();
      Assert (index == triangulation.n_active_cells(), ExcInternalError());

				       // flag all other cells for coarsening
				       // (this should ensure that at least
				       // some of them will actually be
				       // coarsened)
      index=0;
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = triangulation.begin_active();
	   cell != triangulation.end(); ++cell, ++index)
	if (!flags[index])
	  cell->set_coarsen_flag();

      triangulation.execute_coarsening_and_refinement ();
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

  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("crash_03").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog << "OK" << std::endl;
    }
  else
    {
      test<2>();
    }

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
