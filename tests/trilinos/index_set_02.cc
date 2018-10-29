// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



// Test conversion from epetra map back to indexset, this was broken for
// overlapping IndexSets

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"


void
test()
{
  IndexSet set_my(100);
  IndexSet set_ghost(100);

  unsigned int n_proc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  unsigned int myid   = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (myid == 0)
    {
      set_my.add_range(0, 50);
      set_ghost.add_range(0, 50);
      set_ghost.add_range(55, 60);
    }
  else if (myid == 1)
    {
      set_my.add_range(50, 100);
      set_ghost.add_range(45, 100);
    }
  else
    {}

  auto check = [&](IndexSet &idxset) {
    deallog << "IndexSet before size=" << idxset.size() << " values: ";
    idxset.print(deallog);
    IndexSet back(idxset.make_trilinos_map(MPI_COMM_WORLD, true));
    deallog << "IndexSet after size=" << back.size() << " values: ";
    back.print(deallog);
  };

  deallog << "without overlap:" << std::endl;
  check(set_my);

  deallog << "with overlap:" << std::endl;
  check(set_ghost);

  TrilinosWrappers::MPI::Vector v;

  v.reinit(set_my, set_ghost, MPI_COMM_WORLD);
  IndexSet from_partitioner(v.vector_partitioner());
  deallog << "vec size: " << v.size()
          << " from_partitioner: " << from_partitioner.size() << std::endl;

  from_partitioner.print(deallog);
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log;
  test();
}
