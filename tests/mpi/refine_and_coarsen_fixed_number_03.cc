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


// Like _01, but make sure we don't coarsen if we don't want to

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>


#include <fstream>


void test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<2> tr(MPI_COMM_WORLD);

  std::vector<unsigned int> sub(2);
  sub[0] = 5*Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);
  sub[1] = 1;
  GridGenerator::subdivided_hyper_rectangle(static_cast<Triangulation<2>&>(tr),
					    sub, Point<2>(0,0), Point<2>(1,1));
  tr.refine_global (1);

  Vector<float> indicators (tr.dealii::Triangulation<2>::n_active_cells());
  {
    unsigned int cell_index = 0;
    unsigned int my_cell_index = 0;
    for (Triangulation<2>::active_cell_iterator
	   cell = tr.begin_active(); cell != tr.end(); ++cell, ++cell_index)
      if (cell->subdomain_id() == myid)
	{
	  ++my_cell_index;
	  indicators(cell_index) = my_cell_index;
	}
  }


  parallel::distributed::GridRefinement
    ::refine_and_coarsen_fixed_number (tr, indicators, 0.2, 0);

				   // now count number of cells
				   // flagged for refinement and
				   // coarsening. we have to
				   // accumulate over all processors
  unsigned int my_refined   = 0,
	       my_coarsened = 0;
  for (Triangulation<2>::active_cell_iterator
	 cell = tr.begin_active(); cell != tr.end(); ++cell)
    if (cell->refine_flag_set())
      ++my_refined;
    else if (cell->coarsen_flag_set())
      ++my_coarsened;

  unsigned int n_refined   = 0,
	       n_coarsened = 0;
  MPI_Reduce (&my_refined, &n_refined, 1, MPI_UNSIGNED, MPI_SUM, 0,
	      MPI_COMM_WORLD);
  MPI_Reduce (&my_coarsened, &n_coarsened, 1, MPI_UNSIGNED, MPI_SUM, 0,
	      MPI_COMM_WORLD);

				   // make sure we have indeed flagged
				   // exactly 20% of cells for
				   // refinement, and none for
				   // coarsening
  if (myid == 0)
    {
      deallog << "total active cells = "
		<< tr.n_global_active_cells() << std::endl;
      deallog << "n_refined = " << n_refined << std::endl;
      deallog << "n_coarsened = " << n_coarsened << std::endl;
    }

  tr.execute_coarsening_and_refinement ();
  if (myid == 0)
    deallog << "total active cells = "
	    << tr.n_global_active_cells() << std::endl;
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
      std::ofstream logfile(output_file_for_mpi("refine_and_coarsen_fixed_number_03").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();


#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
