//---------------------------------------------------------------------------
//    $Id: 3d_refinement_01.cc 29075 2013-03-27 15:07:57Z heister $
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test that we abort if you refine further than p4est::max_level

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>


template<int dim>
void test(std::ostream& /*out*/)
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  for (unsigned int i=0;i<19;++i)
    {
      deallog << "cells: " << tr.n_active_cells() << " level:" << tr.n_levels() << std::endl;

      typename parallel::distributed::Triangulation<dim>::cell_iterator it;
      it=tr.begin_active();
      while (it->level()<static_cast<int>(i))
	++it;

      
      it->set_refine_flag();
      try
	{
	  tr.execute_coarsening_and_refinement ();
	}
      catch (ExceptionBase &e)
	{
	  deallog << e.get_exc_name() << std::endl;
	}
    }
  
}


int main(int argc, char *argv[])
{
  deal_II_exceptions::disable_abort_on_exception();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  std::ofstream logfile("3d_refinement_09/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();
}
