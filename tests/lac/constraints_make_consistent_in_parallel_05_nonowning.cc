// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test AffineConstraints<double>::make_consistent_in_parallel for a
// case where each process owns one DoF and on each process p>0 we
// know about one constraint that constrains X_p=0.5*X_{p-1}. After making
// things consistent in parallel, this needs to all reduce to X_p=factor*X1.
//
// Compared to the plain _05 test, this variation has the same
// constraints except that the constraint is not added by the process
// that *owns* the DoF, but by another process that only has this DoF
// among its locally relevant ones.

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

  // Process N only owns DoF N:
  const unsigned int n_indices = n_procs;
  IndexSet           owned_elements(n_indices);
  owned_elements.add_index(my_proc);

  // Process N has N, N+1 and N-1 as its locally relevant DoFs:
  IndexSet relevant_elements(n_indices);
  relevant_elements.add_index(my_proc);
  if (my_proc > 0)
    relevant_elements.add_index(my_proc - 1);
  if (my_proc != n_procs - 1)
    relevant_elements.add_index(my_proc + 1);
  AffineConstraints<double> constraints(owned_elements, relevant_elements);

  // Let each process add a constraint for a DoF that it doesn't own:
  if (my_proc != n_procs - 1)
    constraints.add_constraint(my_proc + 1, {{my_proc, 0.5}});

  deallog << "Constraints as created:" << std::endl;
  constraints.print(deallog.get_file_stream());

  constraints.make_consistent_in_parallel(owned_elements,
                                          relevant_elements,
                                          MPI_COMM_WORLD);

  deallog << "Constraints as made consistent:" << std::endl;
  constraints.print(deallog.get_file_stream());

  return 0;
}
