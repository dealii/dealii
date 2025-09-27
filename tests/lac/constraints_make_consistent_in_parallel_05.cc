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



// test AffineConstraints<double>::make_consistent_in_parallel for a
// case where each process owns one DoF and on each process p>0 we
// know about one constraint that constrains X_p=0.5*X_{p-1}. After making
// things consistent in parallel, this needs to all reduce to X_p=factor*X1.

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

  const unsigned int n_indices = n_procs;
  IndexSet           owned_elements(n_indices);
  owned_elements.add_index(my_proc);
  IndexSet relevant_elements(n_indices);
  relevant_elements.add_index(my_proc);
  if (my_proc > 0)
    relevant_elements.add_index(my_proc - 1);
  AffineConstraints<double> constraints(owned_elements, relevant_elements);

  if (my_proc > 0)
    constraints.add_constraint(my_proc, {{my_proc - 1, 0.5}});

  deallog << "Constraints as created:" << std::endl;
  constraints.print(deallog.get_file_stream());

  constraints.make_consistent_in_parallel(owned_elements,
                                          relevant_elements,
                                          MPI_COMM_WORLD);

  deallog << "Constraints as made consistent:" << std::endl;
  constraints.print(deallog.get_file_stream());

  return 0;
}
