// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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

  // Create a copy of v_one using the internal Vec
  PETScWrappers::MPI::Vector v_one_from_vec(v_one.petsc_vector());
  Assert(v_one.size() == v_one_from_vec.size(), ExcInternalError());
  Assert(v_one.locally_owned_size() == v_one_from_vec.locally_owned_size(),
         ExcInternalError());
  Assert(v_one.has_ghost_elements() == v_one_from_vec.has_ghost_elements(),
         ExcInternalError());
  Assert(v_one.ghost_elements() == v_one_from_vec.ghost_elements(),
         ExcInternalError());

  // Test swap
  IndexSet local_active_2(numproc * 4);
  local_active_2.add_range(myid * 4, myid * 4 + 4);
  PETScWrappers::MPI::Vector vs1(local_active, MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector vs2(local_active_2,
                                 local_relevant,
                                 MPI_COMM_WORLD);
  vs1.swap(vs2);
  Assert(vs1.size() == numproc * 4, ExcInternalError());
  Assert(vs2.size() == numproc * 2, ExcInternalError());
  Assert(vs1.locally_owned_size() == 4, ExcInternalError());
  Assert(vs2.locally_owned_size() == 2, ExcInternalError());
  Assert(vs1.has_ghost_elements() == true, ExcInternalError());
  Assert(vs2.has_ghost_elements() == false, ExcInternalError());
  Assert(vs1.ghost_elements() == v_one.ghost_elements(), ExcInternalError());

  // Now BlockVectors
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
  Assert(vb2.petsc_vector() == v.petsc_vector(), ExcInternalError());
  for (unsigned int bl = 0; bl < 2; ++bl)
    {
      Assert(vb2.block(bl).size() == v.block(bl).size(), ExcInternalError());
      Assert(vb2.block(bl).petsc_vector() == v.block(bl).petsc_vector(),
             ExcInternalError());
    }

  // Create new block vector from an array of PETSc vectors
  std::array<Vec, 2> arrayVecs = {
    {vb.block(0).petsc_vector(), vb.block(1).petsc_vector()}};
  PETScWrappers::MPI::BlockVector vb3(arrayVecs);
  Assert(vb3.n_blocks() == vb.n_blocks(), ExcInternalError());
  Assert(vb3.size() == vb.size(), ExcInternalError());
  for (unsigned int bl = 0; bl < 2; ++bl)
    {
      Assert(vb3.block(bl).size() == vb.block(bl).size(), ExcInternalError());
      Assert(vb3.block(bl).petsc_vector() == vb.block(bl).petsc_vector(),
             ExcInternalError());
    }


  // Test swap
  auto old_v_vb2 = vb2.petsc_vector();
  auto old_v_vb  = vb.petsc_vector();
  vb.swap(vb2);
  Assert(vb.petsc_vector() == old_v_vb2, ExcInternalError());
  Assert(vb2.petsc_vector() == old_v_vb, ExcInternalError());

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
