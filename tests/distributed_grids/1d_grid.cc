// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2017 by the deal.II authors
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



// Test that we can compile p::d::Tria<1>

#include "../tests.h"
#include <deal.II/distributed/tria.h>



int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  std::ofstream logfile("output");
  deallog.attach(logfile);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  deallog << "Ok" << std::endl;

}
