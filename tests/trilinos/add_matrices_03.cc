// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2014 by the deal.II authors
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



// Test SparseMatrix::add(factor, SparseMatrix) based on matrices of the same
// sparsity pattern

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <fstream>
#include <iostream>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;


  // each processor owns 3 indices
  IndexSet local_owned(numproc*3);
  local_owned.add_range(myid*3,myid*3+3);

  // Create sparsity patterns
  TrilinosWrappers::SparsityPattern sp(local_owned, MPI_COMM_WORLD);

  for (unsigned int i=myid*3; i<myid*3+3; ++i)
    for (unsigned int j=0; j<local_owned.size(); ++j)
      if ((i+j) % 2 == 1)
        {
          sp.add (i,j);
        }

  sp.compress ();

  // create matrices by adding some elements into the respective positions
  TrilinosWrappers::SparseMatrix m1(sp), m2(sp);
  for (unsigned int i=myid*3; i<myid*3+3; ++i)
    for (unsigned int j=0; j<local_owned.size(); ++j)
      if ((i+j) % 2 == 1)
        {
          m1.add (i,j, i+j);
          if (j%2 == 0)
            m2.add(i,j, i+2*j+1);
        }
  m1.compress(VectorOperation::add);
  m2.compress(VectorOperation::add);

  m1.add(2, m2);

  // Check for correctness of entries (all floating point comparisons should
  // be exact)
  for (unsigned int i=myid*3; i<myid*3+3; ++i)
    for (unsigned int j=0; j<local_owned.size(); ++j)
      if ((i+j) % 2 == 1 && j%2 == 0)
        {
          Assert(m1.el(i,j) == (double)i+j+2*i+4*j+2, ExcInternalError());
        }
      else if ((i+j) % 2 == 1)
        {
          Assert(m1.el(i,j) == (double)i+j, ExcInternalError());
        }
      else
        {
          Assert(m1.el(i,j) == 0., ExcInternalError());
        }

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  try
    {
      test ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
