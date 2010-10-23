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


// like the second part of 3d_refinement_06 but use only four coarse
// grid cells. we want to make sure that the 2:1 relationship holds
// across an edge


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
  {
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

    std::vector<unsigned int> subdivisions(3, 2);
    subdivisions[2] = 1;
    GridGenerator::subdivided_hyper_rectangle(tr, subdivisions, Point<3>(0,0,0), Point<3>(2,2,1));
    tr.begin_active()->set_refine_flag();
    tr.execute_coarsening_and_refinement ();
    for (unsigned int c=0; c<8; ++c)
      tr.begin(0)->child(c)->set_refine_flag();
    tr.execute_coarsening_and_refinement ();

//    write_vtk (tr, "3d_refinement_07", "2");
    deallog << "cells test2: " << tr.n_active_cells() << std::endl;

    Assert (tr.n_active_cells() == 88, ExcInternalError());
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

  std::ofstream logfile("3d_refinement_07/output");
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
