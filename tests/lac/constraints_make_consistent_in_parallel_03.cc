// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
      constraints.reinit(owned_elements, relevant_elements);
    }
  else if (my_proc == 1)
    {
      owned_elements.add_range(1, 3);
      relevant_elements.add_range(0, n_indices);
      constraints.reinit(owned_elements, relevant_elements);
      constraints.add_constraint(0, {{1, 0.5}, {2, 0.5}});
      constraints.add_constraint(2, {{3, 1.5}});
    }
  else if (my_proc == 2)
    {
      owned_elements.add_range(3, n_indices);
      relevant_elements.add_range(0, n_indices);
      constraints.reinit(owned_elements, relevant_elements);
      constraints.add_constraint(0, {{1, 0.5}, {3, 0.75}});
    }

  constraints.make_consistent_in_parallel(owned_elements,
                                          relevant_elements,
                                          MPI_COMM_WORLD);
  constraints.close();
  constraints.print(deallog.get_file_stream());

  return 0;
}
