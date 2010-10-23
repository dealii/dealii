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


// Like coarse_grid_x_01, but instead of checking that what we copy into
// the p4est data structure is correct, check that what we get back is
// correct

#include "../tests.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <grid/tria.h>
#include <distributed/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>

#include <fstream>


template<int dim>
void test(std::ostream& /*out*/)
{
  if (true)
    {
      deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);

      deallog << "Triangulation copied back from p4est has "
	      << tr.n_active_cells()
	      << " active cells" << std::endl;
      GridOut().write_gnuplot (tr, deallog.get_file_stream());
    }


  if (true)
    {
      deallog << "hyper_ball" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_ball(tr, Point<dim>(), 3.);

      deallog << "Triangulation copied back from p4est has "
	      << tr.n_active_cells()
	      << " active cells" << std::endl;
      GridOut().write_gnuplot (tr, deallog.get_file_stream());
    }

  if (true)
    {
      deallog << "half_hyper_ball" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::half_hyper_ball(tr, Point<dim>(), 3.);

      deallog << "Triangulation copied back from p4est has "
	      << tr.n_active_cells()
	      << " active cells" << std::endl;
      GridOut().write_gnuplot (tr, deallog.get_file_stream());
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

  std::ofstream logfile("2d_coarse_grid_x_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
