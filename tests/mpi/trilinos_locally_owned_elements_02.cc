// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2016 by the deal.II authors
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



// test TrilinosVector::locally_owned_elements


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/trilinos_vector.h>

#include <fstream>
#include <sstream>


void test ()
{
  const int n_proc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  //All processes should own 10 entries
  const int entries_per_process = 10;

  IndexSet locally_owned(entries_per_process*n_proc);
  const int begin_index = my_id*entries_per_process;
  const int end_index = (my_id+1)*entries_per_process;
  locally_owned.add_range(begin_index, end_index);

  IndexSet locally_relevant(entries_per_process*n_proc);
  const int local_begin = std::max(0, begin_index-entries_per_process/2);
  const int local_end = entries_per_process*n_proc;
  locally_relevant.add_range (local_begin, local_end);

  TrilinosWrappers::MPI::Vector ghosted, distributed;
  distributed.reinit(locally_owned, MPI_COMM_WORLD);
  ghosted.reinit (locally_owned, locally_relevant, MPI_COMM_WORLD);

  IndexSet locally_owned_elements_distributed = distributed.locally_owned_elements();
  IndexSet locally_owned_elements_ghosted = ghosted.locally_owned_elements();

  const types::global_dof_index local_range_begin_ghosted = ghosted.local_range().first;
  const types::global_dof_index local_range_end_ghosted = ghosted.local_range().second;

  const types::global_dof_index local_range_begin_distributed = distributed.local_range().first;
  const types::global_dof_index local_range_end_distributed = distributed.local_range().second;

  deallog << "locally_owned_elements_distributed: ";
  locally_owned_elements_distributed.print(deallog);
  deallog << "locally_owned_elements_ghosted: ";
  locally_owned_elements_ghosted.print(deallog);
  deallog << "local_range_begin_ghosted: "
          << local_range_begin_ghosted << std::endl;
  deallog << "local_range_end_ghosted: "
          << local_range_end_ghosted << std::endl;
  deallog << "local_range_begin_distributed: "
          << local_range_begin_distributed << std::endl;
  deallog << "local_range_end_distributed: "
          << local_range_end_distributed << std::endl;

  AssertThrow (locally_owned_elements_distributed == locally_owned_elements_ghosted,
               ExcInternalError());

  deallog << "OK" << std::endl;
}



int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  MPILogInitAll log;

  test();
}
