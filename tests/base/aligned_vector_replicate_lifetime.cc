// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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
