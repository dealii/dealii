// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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
  const char * p = strrchr(filename,'/');
  deallog.push (p);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridIn<dim> gi;
  gi.attach_triangulation (tr);
  std::ifstream in (filename);
  gi.read_xda (in);

  write_vtk(tr, "1");

  deallog.pop ();
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("3d");

  test<3> (SOURCE_DIR "/../deal.II/grid_in_3d/1.in");
  test<3> (SOURCE_DIR "/../deal.II/grid_in_3d/2.in");
  test<3> (SOURCE_DIR "/../deal.II/grid_in_3d/3.in");
  test<3> (SOURCE_DIR "/../deal.II/grid_in_3d/4.in");

  test<3> (SOURCE_DIR "/../deal.II/grid_in_3d/evil_0.in");
  test<3> (SOURCE_DIR "/../deal.II/grid_in_3d/evil_1.in");
  test<3> (SOURCE_DIR "/../deal.II/grid_in_3d/evil_2.in");
  test<3> (SOURCE_DIR "/../deal.II/grid_in_3d/evil_3.in");
  test<3> (SOURCE_DIR "/../deal.II/grid_in_3d/evil_4.in");

  deallog.pop();


}
