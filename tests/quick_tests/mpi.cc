// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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

// test that MPI is working correctly. Note that this test expects to
// be executed with exactly two threads.

#include <deal.II/grid/tria.h>
#include <stdio.h>
#include <sched.h>
#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[] )
{
  MPI_Init( &argc, &argv ); 

  int myrank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  std::cout << " Hi from " << myrank << "/" << nproc << std::endl;
  
  if (nproc != 2)
    {
      std::cerr << "ERROR: process does not see nproc=2!" << std::endl;
      return -1;
    }
  
  MPI_Barrier(MPI_COMM_WORLD);

  int err;
  int value = myrank;
  
  if (myrank==1)
    err = MPI_Send(&value, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
  else if (myrank==0)
    err = MPI_Recv(&value, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  if (myrank==0 && value!=1)
    {
      std::cerr << "ERROR: MPI_Send/Recv did not work!" << std::endl;
      return -1;
    }

  value = 1;
  int output = 0;
  
  MPI_Allreduce(&value, &output, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (output != nproc)
    {
      std::cerr << "ERROR: MPI_Allreduce doesn't seem to work!" << std::endl;
      return -1;
    }
    
  // we need this, otherwise gcc will not link against deal.II
  dealii::Triangulation<2> test;

  MPI_Finalize();
}
