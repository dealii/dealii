// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
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

// test ghosted TrilinosWrappers::MPI::Vector -> RWV

#include "../tests.h"
#include <deal.II/base/index_set.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>

#include <vector>

void
test()
{
  IndexSet is(8);
  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank==0)
    is.add_range(0,4);
  if (rank==1)
    is.add_range(4,8);
  is.compress();

  deallog << "IS: ";
  is.print(deallog);

  IndexSet is_ghosted(8);
  is_ghosted.add_range(2,6);
  deallog << "IS_ghosted: ";
  is_ghosted.print(deallog);

  TrilinosWrappers::MPI::Vector tril_vector(is);
  TrilinosWrappers::MPI::Vector tril_vector_ghosted;
  tril_vector_ghosted.reinit(is, is_ghosted, MPI_COMM_WORLD);
  for (unsigned int i=0; i<8; ++i)
    tril_vector[i] = i;

  tril_vector.compress(VectorOperation::insert);
  deallog << "trilinos vec:" << std::endl;
  tril_vector.print(deallog.get_file_stream());

  tril_vector_ghosted = tril_vector;

  deallog << "trilinos vec ghosted:" << std::endl;
  tril_vector_ghosted.print(deallog.get_file_stream());


  IndexSet readwrite_is (tril_vector_ghosted.vector_partitioner());
  deallog << "ghosted IS: ";
  readwrite_is.print(deallog);

  deallog << "tril_vector_ghosted.owned_elements() ";
  tril_vector_ghosted.locally_owned_elements().print(deallog);

  {
    deallog << "import of ghosted vector should fail:" << std::endl;
    LinearAlgebra::ReadWriteVector<double> readwrite(readwrite_is);
    try
      {
        readwrite.import(tril_vector_ghosted, VectorOperation::insert);
        deallog << "RWVector contents:" << std::endl;
        readwrite.print(deallog.get_file_stream());
      }
    catch (ExceptionBase &exc)
      {
        deallog << exc.get_exc_name() << std::endl;
      }
  }
  {
    deallog << "import of distributed vector should work:" << std::endl;
    LinearAlgebra::ReadWriteVector<double> readwrite(readwrite_is);
    readwrite.import(tril_vector, VectorOperation::insert);
    deallog << "RWVector contents:" << std::endl;
    readwrite.print(deallog.get_file_stream());
  }

  deallog << "OK" <<std::endl;
}

int
main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  MPILogInitAll log;

  test();
}
