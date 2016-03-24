// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2016 by the deal.II authors
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


// check SparseMatrix::add(other, factor)

#include "../tests.h"
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  const unsigned int s = 10;
  PETScWrappers::SparseMatrix m(s,s,s);
  for (unsigned int k=0; k<m.m(); ++k)
    for (unsigned int l=0; l<=k; ++l)
      m.set (k,l, k+2*l);

  m.compress (VectorOperation::insert);

  PETScWrappers::SparseMatrix m2(s,s,s);
  m2.set(0,1,5.0);
  for (unsigned int k=0; k<m2.m(); ++k)
    m2.set(k,k,PetscScalar(1.0+k,-1.0-k));
  m2.compress(VectorOperation::insert);

  // we now print the matrix m, it is all real. Then print after
  // adding m2 which is complex, and then subtract m2 again to get the
  // original matrix back.

  deallog << "before: " << m(0,1) << std::endl;
  for (unsigned int k=0; k<s; ++k)
    deallog << m(k,k) << " ";
  deallog << std::endl;

  m.add(m2,1.0);

  deallog << "after: " << m(0,1) << std::endl;
  for (unsigned int k=0; k<s; ++k)
    deallog << m(k,k) << " ";
  deallog << std::endl;

  m.add(m2,-1.0);

  deallog << "back to original: " << m(0,1) << std::endl;
  for (unsigned int k=0; k<s; ++k)
    deallog << m(k,k) << " ";
  deallog << std::endl;

  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        test ();
      }

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
