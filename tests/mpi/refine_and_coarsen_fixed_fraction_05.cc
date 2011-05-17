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


// like _05 but with indicators that span several orders of magnitude

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>


#include <fstream>


void test()
{
  const unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<2> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(static_cast<Triangulation<2>&>(tr));
  tr.refine_global (3);

  Vector<float> indicators (tr.dealii::Triangulation<2>::n_active_cells());
  float min_indicator = tr.dealii::Triangulation<2>::n_active_cells(),
	max_indicator = 0;
  {
    unsigned int cell_index = 0;
    unsigned int my_cell_index = 0;
    for (Triangulation<2>::active_cell_iterator
	   cell = tr.begin_active(); cell != tr.end(); ++cell, ++cell_index)
      if (cell->subdomain_id() == myid)
	{
	  ++my_cell_index;
	  indicators(cell_index) = std::pow(2.,(int)my_cell_index/4);
	  min_indicator = std::min (min_indicator, indicators(cell_index));
	  max_indicator = std::max (max_indicator, indicators(cell_index));
	}
  }

				   // use one strategy to compute
				   // thresholds and obtain those
				   // thresholds
  parallel::distributed::GridRefinement
    ::refine_and_coarsen_fixed_fraction (tr, indicators, 2./3, 1./6);
  {
    float coarsen_indicator = min_indicator-1,
	  refine_indicator  = max_indicator+1;
    unsigned int cell_index=0;
    for (Triangulation<2>::active_cell_iterator
	   cell = tr.begin_active(); cell != tr.end(); ++cell, ++cell_index)
      if (cell->refine_flag_set())
	refine_indicator = std::min (refine_indicator,
				     indicators(cell_index));
      else if (cell->coarsen_flag_set())
	coarsen_indicator = std::max (coarsen_indicator,
				      indicators(cell_index));
    if (myid == 0)
      {
	deallog << "thresholds = " << refine_indicator << ' '
		<< coarsen_indicator << std::endl;
      }

    for (Triangulation<2>::active_cell_iterator
	   cell = tr.begin_active(); cell != tr.end(); ++cell)
      {
	cell->clear_coarsen_flag();
	cell->clear_refine_flag();
      }
  }

				   // now use the second strategy to
				   // compute thresholds and obtain
				   // those thresholds. note that this
				   // only works because we are
				   // working on only a single
				   // processor
  dealii::GridRefinement
    ::refine_and_coarsen_fixed_fraction (tr, indicators, 2./3, 1./6);
  {
    float coarsen_indicator = min_indicator-1,
	  refine_indicator  = max_indicator+1;
    unsigned int cell_index=0;
    for (Triangulation<2>::active_cell_iterator
	   cell = tr.begin_active(); cell != tr.end(); ++cell, ++cell_index)
      if (cell->refine_flag_set())
	refine_indicator = std::min (refine_indicator,
				     indicators(cell_index));
      else if (cell->coarsen_flag_set())
	coarsen_indicator = std::max (coarsen_indicator,
				      indicators(cell_index));
    if (myid == 0)
      {
	deallog << "thresholds = " << refine_indicator << ' '
		<< coarsen_indicator << std::endl;
      }

    for (Triangulation<2>::active_cell_iterator
	   cell = tr.begin_active(); cell != tr.end(); ++cell)
      {
	cell->clear_coarsen_flag();
	cell->clear_refine_flag();
      }
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


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("refine_and_coarsen_fixed_fraction_05").c_str());
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
