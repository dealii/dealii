// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2001 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// check GridTools::distort_random


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>



template <int dim>
void test1 (const bool keep_boundary)
{
  const unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> tria (MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  GridTools::distort_random (0.1, tria, keep_boundary);

  deallog << "dim=" << dim << ", keep_boundary=" << keep_boundary << std::endl;
  std::string filename;
  if (keep_boundary)
    filename = "keep_true-";
  else
    filename = "keep_false-";
  filename += Utilities::int_to_string(dim);

  std::ofstream logfile
    ((filename+ "-" + Utilities::int_to_string(my_id,2)).c_str());

  GridOut().write_gnuplot (tria, logfile);
  MPI_Barrier(MPI_COMM_WORLD);

  if (my_id==0)
    for (unsigned int i=0;
         i<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); ++i)
    {
      deallog << "Process " << i << ":" << std::endl;
      cat_file((filename + "-" + Utilities::int_to_string(i,2)).c_str());
    }

  deallog << "OK" << std::endl;
}



int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, numbers::invalid_unsigned_int);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
  {
    std::ofstream logfile("output");
    deallog << std::setprecision(4);
    logfile << std::setprecision(4);
    deallog.attach(logfile);
    deallog.depth_console(0);
    deallog.threshold_double(1.e-10);

    test1<2> (true);
    test1<2> (false);
    test1<3> (true);
    test1<3> (false);
  }
  else
  {
    test1<2> (true);
    test1<2> (false);
    test1<3> (true);
    test1<3> (false);
  }

  return 0;
}
