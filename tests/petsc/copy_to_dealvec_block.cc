// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2022 by the deal.II authors
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



// Test
// LinearAlgebra::distributed::Vector::operator=(PETScWrappers::MPI::BlockVector&)
// and PETScWrappers::MPI::BlockVector interaction with PETSc Vecs
#include <deal.II/base/index_set.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  // each processor owns 2 indices and all
  // are ghosting Element 1 (the second)

  IndexSet local_active(numproc * 2);
  local_active.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(numproc * 2);
  local_relevant.add_range(1, 2);

  PETScWrappers::MPI::Vector vb_one(local_active, MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector v_one(local_active,
                                   local_relevant,
                                   MPI_COMM_WORLD);

  LinearAlgebra::distributed::Vector<double> copied_one(local_active,
                                                        local_relevant,
                                                        MPI_COMM_WORLD);

  // set local values
  vb_one(myid * 2)     = myid * 2.0;
  vb_one(myid * 2 + 1) = myid * 2.0 + 1.0;

  vb_one.compress(VectorOperation::insert);
  vb_one *= 2.0;
  v_one = vb_one;

  PETScWrappers::MPI::BlockVector vb, v;
  vb.reinit(2);
  v.reinit(2);
  LinearAlgebra::distributed::BlockVector<double> copied(2);
  for (unsigned int bl = 0; bl < 2; ++bl)
    {
      vb.block(bl)     = vb_one;
      v.block(bl)      = v_one;
      copied.block(bl) = copied_one;
    }
  vb.collect_sizes();
  v.collect_sizes();
  copied.collect_sizes();

  copied = vb;

  // check local values
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      deallog << myid * 2 << ':' << copied(myid * 2) << std::endl;
      deallog << myid * 2 + 1 << ':' << copied(myid * 2 + 1) << std::endl;
    }

  for (unsigned int bl = 0; bl < 2; ++bl)
    {
      Assert(copied.block(bl)(myid * 2) == myid * 4.0, ExcInternalError());
      Assert(copied.block(bl)(myid * 2 + 1) == myid * 4.0 + 2.0,
             ExcInternalError());
    }

  copied = v;

  // check ghost values
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "ghost: " << copied(1) << std::endl;
  Assert(copied(1) == 2.0, ExcInternalError());
  Assert(copied.block(1)(1) == 2.0, ExcInternalError());

  // check local values
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      deallog << myid * 2 << ':' << copied(myid * 2) << std::endl;
      deallog << myid * 2 + 1 << ':' << copied(myid * 2 + 1) << std::endl;
    }

  for (unsigned int bl = 0; bl < 2; ++bl)
    {
      Assert(copied.block(bl)(myid * 2) == myid * 4.0, ExcInternalError());
      Assert(copied.block(bl)(myid * 2 + 1) == myid * 4.0 + 2.0,
             ExcInternalError());
    }

  // Create new block vector from a PETSc VECNEST
  PETScWrappers::MPI::BlockVector vb2(v.petsc_vector());
  Assert(vb2.n_blocks() == v.n_blocks(), ExcInternalError());
  Assert(vb2.size() == v.size(), ExcInternalError());
  for (unsigned int bl = 0; bl < 2; ++bl)
    {
      Assert(vb2.block(bl).size() == v.block(bl).size(), ExcInternalError());
      Assert(vb2.block(bl).petsc_vector() == v.block(bl).petsc_vector(),
             ExcInternalError());
    }

  // done
  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
