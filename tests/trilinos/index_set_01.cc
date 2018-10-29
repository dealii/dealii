// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// check that make_trilinos_map always uses the same constructor all the
// processor. If we don't the code will hang like reported on the mailing list
// by Nicola Giuliani.

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <iostream>

#include "../tests.h"



void
test()
{
  IndexSet my_set;
  my_set.set_size(654);
  int this_mpi_process, n_mpi_processes;
  n_mpi_processes  = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  this_mpi_process = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (this_mpi_process == 0)
    {
      my_set.add_range(0, 86);
    }
  else if (this_mpi_process == 1)
    {
      my_set.add_range(86, 400);
      my_set.add_range(529, 654);
    }
  else
    {
      my_set.add_range(400, 529);
    }
  my_set.compress();
  my_set.make_trilinos_map(MPI_COMM_WORLD, false);
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  test();

  deallog << "OK" << std::endl;
}
