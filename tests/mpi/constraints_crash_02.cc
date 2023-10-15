// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
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
  cm.add_line(0);
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
