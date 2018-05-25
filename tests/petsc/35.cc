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



// check PETScWrappers::MPI::Vector::operator+=(Vector)

#include <deal.II/lac/petsc_parallel_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(PETScWrappers::MPI::Vector &v, PETScWrappers::MPI::Vector &w)
{
  // set only certain elements of each
  // vector
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      v(i) = i;
      if (i % 3 == 0)
        w(i) = i + 1.;
    }

  v.compress(VectorOperation::insert);
  w.compress(VectorOperation::insert);

  v += w;

  // make sure we get the expected result
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      if (i % 3 == 0)
        {
          AssertThrow(w(i) == i + 1., ExcInternalError());
          AssertThrow(v(i) == i + i + 1., ExcInternalError());
        }
      else
        {
          AssertThrow(w(i) == 0, ExcInternalError());
          AssertThrow(v(i) == i, ExcInternalError());
        }
    }


  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        IndexSet indices(100);
        indices.add_range(0, 100);
        PETScWrappers::MPI::Vector v(indices, MPI_COMM_WORLD);
        PETScWrappers::MPI::Vector w(indices, MPI_COMM_WORLD);
        test(v, w);
      }
    }
  catch (std::exception &exc)
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
  catch (...)
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
