//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test that anisotropic refinement really doesn't work

#include "../tests.h"
#include "coarse_grid_common.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <grid/tria.h>
#include <distributed/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>

#include <fstream>


template<int dim>
void test(std::ostream& /*out*/)
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.begin_active()->set_refine_flag(RefinementCase<dim>::cut_x);

  try
    {
      tr.execute_coarsening_and_refinement ();
    }
  catch (const ExcMessage &e)
    {
      deallog << "Caught exception:" << std::endl;
      deallog << e.what() << std::endl;
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

  std::ofstream logfile("2d_refinement_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

				   // we want to catch exceptions
				   // instead of aborting the program
  deal_II_exceptions::disable_abort_on_exception();

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();
  deallog.push("2d");
  test<3>(logfile);
  deallog.pop();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
