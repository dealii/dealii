// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// AffineConstraints<double>::distribute() crashes because copy_from() is
// incomplete:

/**
An error occurred in line <2679> of file
</ssd/deal-git1/include/deal.II/lac/affine_constraints.templates.h> in function
    void dealii::AffineConstraints<number>::distribute(VectorType&) const [with
VectorType = dealii::LinearAlgebra::distributed::Vector<double>; number =
double] The violated condition was: needed_elements_for_distribute != IndexSet()
*/

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"



void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  IndexSet owned(numproc);
  owned.add_range(myid, myid + 1);
  IndexSet relevant = complete_index_set(numproc);


  AffineConstraints<double> cm(owned, relevant);
  cm.constrain_dof_to_zero(0);
  cm.close();

  AffineConstraints<double> cm2;
  cm2.copy_from(cm);

  LinearAlgebra::distributed::Vector<double> x(owned, relevant, MPI_COMM_WORLD);

  cm.distribute(x);  // this works
  cm2.distribute(x); // this used to fail
  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test();
  return 0;
}
