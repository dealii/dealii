// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test constructor/reinit of BlockVector with IndexSets

#include <deal.II/base/index_set.h>

#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <class BLOCKVEC>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  std::vector<IndexSet> local_active(2);

  // block 0:
  local_active[0].set_size(numproc);
  local_active[0].add_range(myid, myid + 1);

  local_active[1].set_size(2 * numproc);
  local_active[1].add_range(myid * 2, myid * 2 + 2);

  BLOCKVEC v(local_active, MPI_COMM_WORLD);

  v(myid) = 100.0 + myid;

  v.block(1)(myid * 2)     = myid * 2.0;
  v.block(1)(myid * 2 + 1) = myid * 2.0 + 1.0;

  v.compress(VectorOperation::insert);

  // Assert(!v.has_ghost_elements(), ExcInternalError());

  deallog << "size: " << v.size() << std::endl;
  deallog << "size[0]: " << v.block(0).size() << std::endl;
  deallog << "size[1]: " << v.block(1).size() << std::endl;

  {
    std::ofstream file(
      (std::string("dat.") + Utilities::int_to_string(myid)).c_str());

    file << "**** proc " << myid << std::endl;
    v.print(file);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid == 0)
    {
      for (unsigned int i = 0;
           i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
           ++i)
        {
          cat_file((std::string("dat.") + Utilities::int_to_string(i)).c_str());
        }
    }

  // now test the constructor

  IndexSet id(7);
  if (myid == 0)
    id.add_range(0, 7);
  local_active.push_back(id);

  BLOCKVEC xy(local_active);

  deallog << "size: " << xy.size() << std::endl;
  deallog << "size[0]: " << xy.block(0).size() << std::endl;
  deallog << "size[1]: " << xy.block(1).size() << std::endl;
  deallog << "size[2]: " << xy.block(2).size() << std::endl;

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

      test<PETScWrappers::MPI::BlockVector>();
      // test<TrilinosWrappers::MPI::BlockVector>();
    }
  else
    {
      test<PETScWrappers::MPI::BlockVector>();
      // test<TrilinosWrappers::MPI::BlockVector>();
    }
}
