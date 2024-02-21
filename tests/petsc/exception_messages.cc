// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that we can catch some standard PETSc error codes and print
// appropriate exception messages to the screen.

#include <deal.II/lac/exceptions.h>

#include <petscconf.h>
#include <petscsys.h>

#include <vector>

#include "../tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  std::vector<int> petsc_error_codes;
  petsc_error_codes.push_back(PETSC_ERR_MIN_VALUE);
  petsc_error_codes.push_back(PETSC_ERR_MEM);
  petsc_error_codes.push_back(PETSC_ERR_LIB);
  petsc_error_codes.push_back(PETSC_ERR_ARG_WRONGSTATE);
  petsc_error_codes.push_back(PETSC_ERR_ARG_NOTSAMETYPE);

  for (unsigned int i = 0; i < petsc_error_codes.size(); ++i)
    {
      try
        {
          throw ExcPETScError(petsc_error_codes[i]);
        }
      catch (const ExceptionBase &exc)
        {
          deallog << exc.get_exc_name() << std::endl;
          exc.print_info(deallog.get_file_stream());
        }
    }

  return 0;
}
