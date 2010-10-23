//---------------------------------------------------------------------------
//    $Id$
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


// Test interaction with p4est with a simple grid in 3d. here, we
// check that if we refine a square once, and then one of the children
// once again, that we get 15 cells
//
// at the time of writing this test, the results for this testcase
// were erratic and apparently non-deterministic. the actual cause was
// an uninitialized variable, fixed in revision 16414.

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
#include <grid/grid_in.h>

#include <fstream>
#include <cstdlib>


template<int dim>
void test(std::ostream& /*out*/)
{
  for (unsigned int i=0; i<GeometryInfo<dim>::max_children_per_cell; ++i)
    {
      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);

      deallog << i << ' ' << tr.n_active_cells()
	      << std::endl;

      tr.refine_global (1);

      deallog << i << ' ' << tr.n_active_cells()
	      << std::endl;

      Assert (tr.n_active_cells() == 8,
	      ExcInternalError());

      typename Triangulation<dim>::active_cell_iterator
	cell = tr.begin_active();
      std::advance (cell, i);
      cell->set_refine_flag();
      tr.execute_coarsening_and_refinement ();

      deallog << i << ' ' << tr.n_active_cells()
	      << std::endl;

      Assert (tr.n_active_cells() == 15,
	      ExcInternalError());
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

  std::ofstream logfile("3d_refinement_04/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
