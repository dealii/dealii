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



// verify that VectorBase::print uses the precision parameter correctly and
// restores the previous value of the stream precision

#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        IndexSet indices(5);
        indices.add_range(0, 5);
        PETScWrappers::MPI::Vector v(indices, MPI_COMM_WORLD);
        for (unsigned int i = 0; i < v.size(); ++i)
          v(i) = i * 1.2345678901234567;

        // print with old precision
        deallog << numbers::PI << std::endl;

        // print with prescribed precision
        deallog << "across=false" << std::endl;
        v.print(deallog.get_file_stream(), 10, true, false);

        deallog << "across=true" << std::endl;
        v.print(deallog.get_file_stream(), 10, true, true);

        // print once more. should be the old precision again
        deallog << numbers::PI << std::endl;
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
