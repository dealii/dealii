// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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



// check ConstraintMatrix::set_zero(Vector) for parallel vectors.
// this documents a bug introduced in r29678 for block vectors
// that was fixed in 29940.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>

#include <fstream>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  std::vector<IndexSet> local_active(2);

  // block 0:
  local_active[0].set_size(numproc);
  local_active[0].add_range(myid,myid+1);

  // block 1:
  local_active[1].set_size(2*numproc);
  local_active[1].add_range(myid*2,myid*2+2);

  PETScWrappers::MPI::BlockVector v(local_active, MPI_COMM_WORLD);

  for (unsigned int i=0; i<v.size(); ++i)
    v(i)=1.0+i;
  v.compress(VectorOperation::insert);

  IndexSet local_active_together(3*numproc);
  local_active_together.add_range(myid,myid+1);
  local_active_together.add_range(numproc+myid*2,numproc+myid*2+2);

  ConstraintMatrix cm(local_active_together);
  cm.add_line(numproc + myid*2);
  cm.close();

  deallog << "vector before:" << std::endl;
  v.print(deallog.get_file_stream());

  deallog << "CM:" << std::endl;
  cm.print(deallog.get_file_stream());

  cm.set_zero(v);

  deallog << "vector after:" << std::endl;
  v.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  test();
  return 0;
}
