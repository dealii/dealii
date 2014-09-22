// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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



// test the various VECTOR::extract_subvector_to functions for
// parallel vectors and block vectors

#include "../tests.h"
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>


template <class Vector>
void set (Vector &vector)
{
  for (unsigned int i=0; i<vector.size(); ++i)
    if (vector.locally_owned_elements().is_element(i))
      vector(i) = i;
  vector.compress (VectorOperation::insert);
}


template <class Vector>
void test (Vector &vector)
{
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  // select every other element
  std::vector<typename Vector::size_type> indices;
  for (unsigned int j=0; j<vector.size()/2; ++j)
    indices.push_back (2*j);

  // do the extraction with the function that takes indices, then
  // assert correctness
  std::vector<typename Vector::value_type> values1 (indices.size());
  vector.extract_subvector_to (indices, values1);
  for (unsigned int j=0; j<vector.size()/2; ++j)
    Assert (values1[j] == 2*j, ExcInternalError());

  // do the same with the version of the function that takes iterators
  std::vector<typename Vector::value_type> values2 (indices.size());
  vector.extract_subvector_to (indices.begin(),
			       indices.end(),
			       values2.begin());
  for (unsigned int j=0; j<vector.size()/2; ++j)
    Assert (values2[j] == 2*j, ExcInternalError());

  // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  {
    IndexSet local (10);
    if (myid==0)
      local.add_range(0,5);
    if (myid==1)
      local.add_range(5,10);

    IndexSet dense_local (10);
    dense_local.add_range(0,10);

    {
      deallog.push("deal.II");
      parallel::distributed::Vector<double> w(local, MPI_COMM_WORLD);
      set (w);
      parallel::distributed::Vector<double> v(local, dense_local, MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test (v);
      deallog.pop();
    }

    {
      deallog.push("PETSc");
      PETScWrappers::MPI::Vector w(local, MPI_COMM_WORLD);
      set (w);
      PETScWrappers::MPI::Vector v(local, dense_local, MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test (v);
      deallog.pop();
    }

    {
      deallog.push("Trilinos");
      TrilinosWrappers::MPI::Vector w(local, MPI_COMM_WORLD);
      set (w);
      TrilinosWrappers::MPI::Vector v(local, dense_local, MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test (v);
      deallog.pop();
    }


    std::vector<IndexSet> partitioning;
    {
      IndexSet block1(10);
      if (myid==0)
	block1.add_range(0,7);
      if (myid==1)
	block1.add_range(7,10);

      IndexSet block2(6);
      if (myid==0)
	block2.add_range(0,2);
      if (myid==1)
	block2.add_range(2,6);

      partitioning.push_back(block1);
      partitioning.push_back(block2);
    }

    std::vector<IndexSet> dense_partitioning;
    {
      IndexSet block1(10);
      block1.add_range(0,10);

      IndexSet block2(6);
      block2.add_range(0,6);

      dense_partitioning.push_back(block1);
      dense_partitioning.push_back(block2);
    }


    {
      deallog.push("deal.II");
      parallel::distributed::BlockVector<double> w(partitioning, MPI_COMM_WORLD);
      set (w);
      parallel::distributed::BlockVector<double> v(partitioning, dense_partitioning, MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test (v);
      deallog.pop();
    }

    {
      deallog.push("PETSc");
      PETScWrappers::MPI::BlockVector w(partitioning, MPI_COMM_WORLD);
      set (w);
      PETScWrappers::MPI::BlockVector v(partitioning, dense_partitioning, MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test (v);
      deallog.pop();
    }

    {
      deallog.push("Trilinos");
      TrilinosWrappers::MPI::BlockVector w(partitioning, MPI_COMM_WORLD);
      set (w);
      TrilinosWrappers::MPI::BlockVector v(partitioning, dense_partitioning, MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test (v);
      deallog.pop();
    }

  }

}
