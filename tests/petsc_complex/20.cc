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



// check PETScWrappers::Vector::operator*(Vector) on two vectors that are
// not orthogonal

#include "../tests.h"
#include <deal.II/lac/petsc_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test (PETScWrappers::Vector &v,
           PETScWrappers::Vector &w)
{
  // set only certain elements of each
  // vector, and record the expected scalar
  // product
  PetscScalar product = PetscScalar(0);
  for (unsigned int i=0; i<v.size(); ++i)
    {
      const PetscScalar vi = std::complex<double> (i, (i%2)*2.0);
      v(i) = vi;
      if (i%3 == 0)
        {
          const PetscScalar wi = std::complex<double> (5.0-i,2.5*(i%6));
          w(i) = wi;
          product += PetscConj(vi)*wi;
        }
    }

  v.compress (VectorOperation::insert);
  w.compress (VectorOperation::insert);

  // make sure the scalar product is zero
  AssertThrow (v*w == product, ExcInternalError());

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
        PETScWrappers::Vector v (100);
        PETScWrappers::Vector w (100);
        test (v,w);
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
