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



// check PETScWrappers::SparseMatrix::operator *=

#include "../tests.h"
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <fstream>
#include <iostream>


void test (PETScWrappers::SparseMatrix &m)
{
  // first set a few entries
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        m.set (i,j, std::complex<double> (i*j*.5+.5,i*j*.5));

  m.compress (VectorOperation::insert);

  // then multiply everything by 1.25
  m *= 1.25;

  // and make sure we retrieve the values we expect
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        {
          AssertThrow (m(i,j)    == std::complex<double> ((i*j*.5+.5)*1.25,i*j*.5*1.25), ExcInternalError());
          AssertThrow (m.el(i,j) == std::complex<double> ((i*j*.5+.5)*1.25,i*j*.5*1.25), ExcInternalError());
        }
      else
        {
          AssertThrow (m(i,j)    == std::complex<double> (0.,0.), ExcInternalError());
          AssertThrow (m.el(i,j) == std::complex<double> (0.,0.), ExcInternalError());
        }

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        PETScWrappers::SparseMatrix m(5,5,3);
        test (m);
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


