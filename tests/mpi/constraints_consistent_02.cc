// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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



// test AffineConstraints<double>::make_consistent_in_parallel for case
// where constraints need to be combined between different ranks

#include <deal.II/base/mpi.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  MPILogInitAll                    all;

  const unsigned int my_proc = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  AssertThrow(n_procs == 3, ExcMessage("Only working on 3 ranks"));

  const unsigned int        n_indices = 5;
  IndexSet                  owned_elements(n_indices);
  IndexSet                  relevant_elements(n_indices);
  AffineConstraints<double> constraints;
  if (my_proc == 0)
    {
      owned_elements.add_index(0);
      relevant_elements.add_range(0, 3);
    }
  else if (my_proc == 1)
    {
      owned_elements.add_range(1, 3);
      relevant_elements.add_range(0, n_indices);
      constraints.add_line(0);
      constraints.add_entry(0, 1, 0.5);
      constraints.add_entry(0, 2, 0.5);
      constraints.add_line(2);
      constraints.add_entry(2, 3, 1.5);
    }
  else if (my_proc == 2)
    {
      owned_elements.add_range(3, n_indices);
      relevant_elements.add_range(0, n_indices);
      constraints.add_line(0);
      constraints.add_entry(0, 1, 0.5);
      constraints.add_entry(0, 3, 0.5);
    }

  deallog << "Output 0" << std::endl;
  constraints.print(deallog.get_file_stream());
  constraints.make_consistent_in_parallel(owned_elements,
                                          relevant_elements,
                                          MPI_COMM_WORLD);
  deallog << "Output 1" << std::endl;
  constraints.print(deallog.get_file_stream());
  constraints.close();
  deallog << "Output 2" << std::endl;
  constraints.print(deallog.get_file_stream());

  return 0;
}
