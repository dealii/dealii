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


#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <fstream>
#include <iostream>
#include <vector>

// Check LinearAlgebra::EpetraWrappers::Vector assignement and the operator
// overloading.


void test()
{
  IndexSet parallel_partitioner_1(10);
  IndexSet parallel_partitioner_2(10);
  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank==0)
    {
      parallel_partitioner_1.add_range(0,5);
      parallel_partitioner_2.add_range(0,3);
    }
  else
    {
      parallel_partitioner_1.add_range(5,10);
      parallel_partitioner_2.add_range(3,10);
    }
  parallel_partitioner_1.compress();
  parallel_partitioner_2.compress();
  LinearAlgebra::EpetraWrappers::Vector a;
  LinearAlgebra::EpetraWrappers::Vector b(parallel_partitioner_1, MPI_COMM_WORLD);
  LinearAlgebra::EpetraWrappers::Vector c(b);

  AssertThrow(a.size()==0, ExcMessage("Vector has the wrong size."));
  AssertThrow(b.size()==10, ExcMessage("Vector has the wrong size."));
  AssertThrow(c.size()==10, ExcMessage("Vector has the wrong size."));

  a.reinit(parallel_partitioner_2, MPI_COMM_WORLD);
  AssertThrow(a.size()==10, ExcMessage("Vector has the wrong size."));

  AssertThrow(parallel_partitioner_1==b.locally_owned_elements(),
              ExcMessage("IndexSet has been modified."));
  AssertThrow(parallel_partitioner_2==a.locally_owned_elements(),
              ExcMessage("IndexSet has been modified."));

  IndexSet read_write_index_set(10);
  if (rank==0)
    read_write_index_set.add_range(0,5);
  else
    read_write_index_set.add_range(5,10);
  read_write_index_set.compress();

  LinearAlgebra::ReadWriteVector<double> read_write_1(read_write_index_set);
  LinearAlgebra::ReadWriteVector<double> read_write_2(read_write_index_set);
  LinearAlgebra::ReadWriteVector<double> read_write_3(read_write_index_set);
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        {
          read_write_1[i] = i;
          read_write_2[i] = 5.+i;
        }
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        {
          read_write_1[i] = i;
          read_write_2[i] = 5.+i;
        }
    }

  a = read_write_2;
  b = read_write_1;
  c = read_write_2;

  read_write_3 = a;
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(read_write_2[i]==read_write_3[i],
                    ExcMessage("Vector a has been modified."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(read_write_2[i]==read_write_3[i],
                    ExcMessage("Vector a has been modified."));
    }

  read_write_3 = b;
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(read_write_1[i]==read_write_3[i],
                    ExcMessage("Vector b has been modified."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(read_write_1[i]==read_write_3[i],
                    ExcMessage("Vector b has been modified."));
    }

  read_write_3 = c;
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(read_write_2[i]==read_write_3[i],
                    ExcMessage("Vector c has been modified."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(read_write_2[i]==read_write_3[i],
                    ExcMessage("Vector c has been modified."));
    }


  a *= 2;
  read_write_3 = a;
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(2.*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in operator *=."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(2.*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in operator *=."));
    }

  c /= 2.;
  read_write_3 = c;
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(0.5*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in operator /=."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(0.5*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in operator /=."));
    }

  b += a;
  read_write_3 = b;
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(2.*read_write_2[i]+read_write_1[i]==read_write_3[i],
                    ExcMessage("Problem in operator +=."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(2.*read_write_2[i]+read_write_1[i]==read_write_3[i],
                    ExcMessage("Problem in operator +=."));
    }

  b -= c;
  read_write_3 = b;
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(1.5*read_write_2[i]+read_write_1[i]==read_write_3[i],
                    ExcMessage("Problem in operator -=."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(1.5*read_write_2[i]+read_write_1[i]==read_write_3[i],
                    ExcMessage("Problem in operator -=."));
    }

  b = read_write_1;
  c = read_write_1;
  const double val = b*c;
  AssertThrow(val==285., ExcMessage("Problem in operator *."));
}


int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  deallog << "OK" <<std::endl;

  return 0;
}
