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



// this is equivalent to the petsc_parallel_vector_assign_01 test, except that
// we use operator+= instead of operator=. This is also not exciting...

#include <deal.II/lac/petsc_parallel_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(PETScWrappers::MPI::Vector &v, PETScWrappers::MPI::Vector &w)
{
  // set the first vector
  for (unsigned int k = 0; k < v.size(); ++k)
    v(k) = PetscScalar(k, k * k);

  // add elements by reference
  for (unsigned int k = 0; k < v.size(); ++k)
    w(k) += v(k);

  // check that they're equal
  AssertThrow(v == w, ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        PETScWrappers::MPI::Vector v(MPI_COMM_WORLD, 20, 20);
        PETScWrappers::MPI::Vector w(MPI_COMM_WORLD, 20, 20);
        test(v, w);

        // Some output
        deallog << "Complex vectors: " << std::endl;
        deallog << "v: " << std::flush;
        v.print(logfile, 0, false, true);
        deallog << "w: " << std::flush;
        w.print(logfile, 0, false, true);
        deallog << "OK" << std::endl;
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
