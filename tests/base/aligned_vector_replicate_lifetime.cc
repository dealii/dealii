// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test AlignedVector::replicate_across_communicator().
//
// Check that all objects are correctly created/destroyed again.

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/logstream.h>

#include "../tests.h"


int object_number     = 0;
int objects_destroyed = 0;

class C
{
public:
  C()
  {
    object_number = ::object_number++;
    deallog << "Default constructor. Object number " << object_number
            << std::endl;
  }

  C(const C &c)
  {
    object_number = ::object_number++;
    deallog << "copy constructor from " << c.object_number << ". Object number "
            << object_number << std::endl;
  }

  C(const C &&c)
  {
    object_number = ::object_number++;
    deallog << "move constructor from " << c.object_number << ". Object number "
            << object_number << std::endl;
  }

  ~C()
  {
    deallog << "destructor. Object number " << object_number << std::endl;
    ++objects_destroyed;
  }

  template <typename Archive>
  void
  serialize(Archive &ar, const unsigned int)
  {
    ar &object_number;
  }


private:
  unsigned int object_number;
};



void
test()
{
  const MPI_Comm     communicator = MPI_COMM_WORLD;
  const unsigned int root         = 1;
  Assert(root < Utilities::MPI::n_mpi_processes(communicator),
         ExcInternalError());

  deallog << "*** On process " << Utilities::MPI::this_mpi_process(communicator)
          << ": Starting up..." << std::endl;

  // Create an object of nonzero size and then replicate it.
  AlignedVector<C> avec(Utilities::MPI::this_mpi_process(communicator));

  deallog << "*** On process " << Utilities::MPI::this_mpi_process(communicator)
          << ": Replicating..." << std::endl;

  avec.replicate_across_communicator(communicator, root);

  deallog << "*** On process " << Utilities::MPI::this_mpi_process(communicator)
          << ": Going out of scope" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test();

  deallog << "Objects created: " << object_number << std::endl;
  deallog << "Objects destroyed: " << objects_destroyed << std::endl;
}
