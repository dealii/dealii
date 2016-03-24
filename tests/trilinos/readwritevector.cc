// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

// test operator =

#include "../tests.h"
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>

#include <fstream>
#include <vector>

void test()
{
  IndexSet is(8);
  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank==0)
    is.add_range(0,4);
  if (rank==1)
    is.add_range(4,8);
  is.compress();
  TrilinosWrappers::MPI::Vector tril_vector(is);
  Vector<double> tmp(8);
  for (unsigned int i=0; i<8; ++i)
    tmp[i] = i;
  tril_vector = tmp;

  tril_vector.compress(VectorOperation::insert);

  MPI_Barrier(MPI_COMM_WORLD);

  IndexSet readwrite_is(8);
  if (rank==0)
    {
      readwrite_is.add_range(0,2);
      readwrite_is.add_range(6,8);
    }
  if (rank==1)
    readwrite_is.add_range(2,6);
  readwrite_is.compress();
  LinearAlgebra::ReadWriteVector<double> readwrite(readwrite_is);
  readwrite.import(tril_vector,VectorOperation::insert);
  if (rank==0)
    {
      std::vector<double> comp(4);
      comp[0] = 0.;
      comp[1] = 1.;
      comp[2] = 6.;
      comp[3] = 7.;
      for (unsigned int i=0; i<4; ++i)
        AssertThrow(readwrite.local_element(i)==comp[i],
                    ExcMessage("Element not copied correctly"));
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank==1)
    {
      std::vector<double> comp(4);
      comp[0] = 2.;
      comp[1] = 3.;
      comp[2] = 4.;
      comp[3] = 5.;
      for (unsigned int i=0; i<4; ++i)
        AssertThrow(readwrite.local_element(i)==comp[i],
                    ExcMessage("Element not copied correctly"));
    }

  readwrite.import(tril_vector,VectorOperation::add);

  if (rank==0)
    {
      std::vector<double> comp(4);
      comp[0] = 0.;
      comp[1] = 2.;
      comp[2] = 12.;
      comp[3] = 14.;
      for (unsigned int i=0; i<4; ++i)
        AssertThrow(readwrite.local_element(i)==comp[i],
                    ExcMessage("Element not copied correctly"));
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank==1)
    {
      std::vector<double> comp(4);
      comp[0] = 4.;
      comp[1] = 6.;
      comp[2] = 8.;
      comp[3] = 10.;
      for (unsigned int i=0; i<4; ++i)
        AssertThrow(readwrite.local_element(i)==comp[i],
                    ExcMessage("Element not copied correctly"));
    }

  deallog << "OK" <<std::endl;
}

int main (int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());
  test();
}
