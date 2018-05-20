// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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

// Test the constructor PETScWrappers::VectorBase(const Vec &) that takes an
// existing PETSc vector for complex values.

#include "../tests.h"
#include <deal.II/lac/petsc_parallel_vector.h>
#include <iostream>
#include <vector>

void
test(PETScWrappers::VectorBase& v, PETScWrappers::MPI::Vector& w)
{
  // set the first vector
  for(unsigned int k = 0; k < v.size(); ++k)
    v(k) = std::complex<double>(k, 0.5 * k);

  // copy elements by reference
  for(unsigned int k = 0; k < v.size(); ++k)
    w(k) = v(k);

  // check that they're equal
  AssertThrow(v == w, ExcInternalError());

  deallog << "OK" << std::endl;
}

int
main(int argc, char** argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      Vec                              vpetsc;
      int ierr = VecCreateSeq(PETSC_COMM_SELF, 100, &vpetsc);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
      {
        PETScWrappers::VectorBase  v(vpetsc);
        PETScWrappers::MPI::Vector w(PETSC_COMM_SELF, 100, 100);
        test(v, w);
      }

#if DEAL_II_PETSC_VERSION_LT(3, 2, 0)
      ierr = VecDestroy(vpetsc);
#else
      ierr = VecDestroy(&vpetsc);
#endif

      AssertThrow(ierr == 0, ExcPETScError(ierr));
    }
  catch(std::exception& exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch(...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
