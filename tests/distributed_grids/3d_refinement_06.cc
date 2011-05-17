//---------------------------------------------------------------------------
//    $Id: 3d_refinement_01.cc 17444 2008-10-31 19:35:14Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check if p4est does limit_level_difference_at_vertices in one 3d tree
// and in different trees
// test1 divides the lower-right cell of a square three times
// test2 does the same with a subdivided_hyper_cube
//
// this test is exactly analogous to 2d_refinement_06


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
  {
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

    GridGenerator::hyper_cube(tr);
    tr.begin_active()->set_refine_flag();
    tr.execute_coarsening_and_refinement ();
    tr.begin_active()->set_refine_flag();
    tr.execute_coarsening_and_refinement ();

				     // for better visibility, refine
				     // all children of the sole
				     // level-1 cell once more. this
				     // will introduce more level-1
				     // cells as well.
    for (unsigned int c=0; c<8; ++c)
      tr.begin(1)->child(c)->set_refine_flag();
    tr.execute_coarsening_and_refinement ();

//    write_vtk (tr, "3d_refinement_06", "1");
    deallog << "cells test1: " << tr.n_active_cells() << std::endl;

    Assert (tr.n_active_cells() == 120, ExcInternalError());
  }
  {
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

    GridGenerator::subdivided_hyper_cube(tr, 2);
    tr.begin_active()->set_refine_flag();
    tr.execute_coarsening_and_refinement ();
    for (unsigned int c=0; c<8; ++c)
      tr.begin(0)->child(c)->set_refine_flag();
    tr.execute_coarsening_and_refinement ();

//    write_vtk (tr, "3d_refinement_06", "2");
    deallog << "cells test2: " << tr.n_active_cells() << std::endl;

    Assert (tr.n_active_cells() == 120, ExcInternalError());
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

  std::ofstream logfile("3d_refinement_06/output");
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
