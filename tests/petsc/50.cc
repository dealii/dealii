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

// check PETScWrappers::MPI::Vector::operator = (Vector<T>) with T!=PetscScalar

#include "../tests.h"
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/vector.h>
#include <iostream>
#include <vector>

void
test(PETScWrappers::MPI::Vector& v)
{
  Vector<double> w(v.size());
  Vector<float>  x(v.size());

  for(unsigned int i = 0; i < w.size(); ++i)
    {
      w(i) = i;
      x(i) = i + 1;
    }

  // first copy from w and make sure we get
  // the expected result. then copy from x
  // and do the same. in at least one of the
  // two cases, the template argument to
  // Vector<T> must be different from
  // PetscScalar
  v = w;
  for(unsigned int i = 0; i < v.size(); ++i)
    {
      AssertThrow(w(i) == i, ExcInternalError());
      AssertThrow(v(i) == i, ExcInternalError());
      AssertThrow(x(i) == i + 1, ExcInternalError());
    }

  v = x;
  for(unsigned int i = 0; i < v.size(); ++i)
    {
      AssertThrow(w(i) == i, ExcInternalError());
      AssertThrow(v(i) == i + 1, ExcInternalError());
      AssertThrow(x(i) == i + 1, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}

int
main(int argc, char** argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        IndexSet indices(100);
        indices.add_range(0, 100);
        PETScWrappers::MPI::Vector v(indices, MPI_COMM_WORLD);
        test(v);
      }
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
