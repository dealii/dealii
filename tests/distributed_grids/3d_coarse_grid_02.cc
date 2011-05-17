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


// Test interaction with p4est with a complicated 3d grid read from file.
//
// the files we read here are the ones already used in deal.II/grid_in_3d

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>



template<int dim>
void test(const char *filename)
{
  deallog.push (filename);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridIn<dim> gi;
  gi.attach_triangulation (tr);
  std::ifstream in (filename);
  gi.read_xda (in);

  write_vtk (tr, "3d_coarse_grid_02", "1");

  deallog.pop ();
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  std::ofstream logfile("3d_coarse_grid_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("3d");

  test<3> ("../deal.II/grid_in_3d/1.in");
  test<3> ("../deal.II/grid_in_3d/2.in");
  test<3> ("../deal.II/grid_in_3d/3.in");
  test<3> ("../deal.II/grid_in_3d/4.in");

  test<3> ("../deal.II/grid_in_3d/evil_0.in");
  test<3> ("../deal.II/grid_in_3d/evil_1.in");
  test<3> ("../deal.II/grid_in_3d/evil_2.in");
  test<3> ("../deal.II/grid_in_3d/evil_3.in");
  test<3> ("../deal.II/grid_in_3d/evil_4.in");

  deallog.pop();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
