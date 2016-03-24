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

// Check LinearAlgebra::EpetraWrappers::Vector add and sadd.

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
  LinearAlgebra::EpetraWrappers::Vector a(parallel_partitioner_1, MPI_COMM_WORLD);
  LinearAlgebra::EpetraWrappers::Vector b(parallel_partitioner_1, MPI_COMM_WORLD);
  LinearAlgebra::EpetraWrappers::Vector c(parallel_partitioner_2, MPI_COMM_WORLD);

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

  a.import(read_write_1, VectorOperation::insert);
  b.import(read_write_2, VectorOperation::insert);
  c.import(read_write_2, VectorOperation::insert);

  a.add(1.);
  read_write_3.import(a, VectorOperation::insert);
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(1.+read_write_1[i]==read_write_3[i],
                    ExcMessage("Problem in add(scalar)."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(1.+read_write_1[i]==read_write_3[i],
                    ExcMessage("Problem in add(scalar)."));
    }

  a.add(2.,b);
  read_write_3.import(a, VectorOperation::insert);
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(1.+read_write_1[i]+2.*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in add(scalar,Vector)."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(1.+read_write_1[i]+2.*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in add(scalar,Vector)."));
    }


  LinearAlgebra::EpetraWrappers::Vector d(a);
  a.add(2.,b,3.,d);
  read_write_3.import(a, VectorOperation::insert);
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(4.+4.*read_write_1[i]+10.*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in add(scalar,Vector,scalar,Vector)."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(4.+4.*read_write_1[i]+10.*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in add(scalar,Vector,scalar,Vector)."));
    }


  a.import(read_write_1, VectorOperation::insert);
  a.sadd(3.,2.,c);
  read_write_3.import(a, VectorOperation::insert);
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(3.*read_write_1[i]+2.*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in sadd(scalar,scalar,Vector)."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(3.*read_write_1[i]+2.*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in sadd(scalar,scalar,Vector)."));
    }


  a.import(read_write_1, VectorOperation::insert);
  a.scale(b);
  read_write_3.import(a, VectorOperation::insert);
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(read_write_1[i]*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in scale."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(read_write_1[i]*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in scale."));
    }


  a.equ(2.,c);
  read_write_3.import(a, VectorOperation::insert);
  if (rank==0)
    {
      for (unsigned int i=0; i<5; ++i)
        AssertThrow(2.*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in scale."));
    }
  else
    {
      for (unsigned int i=5; i<10; ++i)
        AssertThrow(2.*read_write_2[i]==read_write_3[i],
                    ExcMessage("Problem in scale."));
    }


  AssertThrow(b.l1_norm()==95., ExcMessage("Problem in l1_norm."));

  const double eps=1e-6;
  AssertThrow(std::fabs(b.l2_norm()-31.3847096)<eps,
              ExcMessage("Problem in l2_norm"));

  AssertThrow(b.linfty_norm()==14., ExcMessage("Problem in linfty_norm."));

  a.import(read_write_1, VectorOperation::insert);
  const double val = a.add_and_dot(2.,a,b);
  AssertThrow(val==1530., ExcMessage("Problem in add_and_dot"));
}


int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  test();

  deallog << "OK" <<std::endl;

  return 0;
}
