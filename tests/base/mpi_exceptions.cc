// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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


// Check that we can catch all MPI-1 error codes and print appropriate
// exception messages to the screen.

#include <deal.II/base/exceptions.h>

#include <mpi.h>

#include <vector>

#include "../tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  std::vector<int> mpi_error_codes;
  mpi_error_codes.push_back(MPI_ERR_BUFFER);
  mpi_error_codes.push_back(MPI_ERR_COUNT);
  mpi_error_codes.push_back(MPI_ERR_TYPE);
  mpi_error_codes.push_back(MPI_ERR_TAG);
  mpi_error_codes.push_back(MPI_ERR_COMM);
  mpi_error_codes.push_back(MPI_ERR_RANK);
  mpi_error_codes.push_back(MPI_ERR_REQUEST);
  mpi_error_codes.push_back(MPI_ERR_ROOT);
  mpi_error_codes.push_back(MPI_ERR_GROUP);
  mpi_error_codes.push_back(MPI_ERR_OP);
  mpi_error_codes.push_back(MPI_ERR_TOPOLOGY);
  mpi_error_codes.push_back(MPI_ERR_UNKNOWN);
  mpi_error_codes.push_back(MPI_ERR_OTHER);
  mpi_error_codes.push_back(MPI_ERR_INTERN);
  mpi_error_codes.push_back(MPI_ERR_LASTCODE);

  for (unsigned int i = 0; i < mpi_error_codes.size(); ++i)
    {
      try
        {
          throw ExcMPI(mpi_error_codes[i]);
        }
      catch (const ExceptionBase &exc)
        {
          deallog << exc.get_exc_name() << std::endl;
          exc.print_info(deallog.get_file_stream());
        }
    }


  return 0;
}
