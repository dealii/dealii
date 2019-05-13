// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// check PETScWrappers::MPI::Vector::operator==(PETScWrappers::MPI::Vector)
// for vectors that are equal

#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(PETScWrappers::MPI::Vector &v, PETScWrappers::MPI::Vector &w)
{
  // set only certain elements of each vector
  for (unsigned int k = 0; k < v.size(); ++k)
    {
      v(k) = PetscScalar(k, k);
      if (k % 3 == 0)
        w(k) = std::complex<double>(k + 1., k + 1);
    }
  v.compress(VectorOperation::insert);
  w.compress(VectorOperation::insert);

  // then copy elements and make sure the vectors are actually equal
  v = w;
  AssertThrow(v == w, ExcInternalError());

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
        PETScWrappers::MPI::Vector v(MPI_COMM_WORLD, 20, 20);
        PETScWrappers::MPI::Vector w(MPI_COMM_WORLD, 20, 20);
        test(v, w);

        // Output
        deallog << "Complex vectors: " << std::endl;
        v.print(deallog.get_file_stream(), 0, false, true);
        w.print(deallog.get_file_stream(), 0, false, true);
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
