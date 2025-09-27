// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check is_serial_vector type trait

#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>

#include "../tests.h"

void
test()
{
  // make sure that is_serial_vector< dealii::PETScWrappers::MPI::Vector > is
  // working
  Assert(is_serial_vector<dealii::PETScWrappers::MPI::Vector>::value == false,
         ExcInternalError());

  deallog << is_serial_vector<dealii::PETScWrappers::MPI::Vector>::value
          << std::endl;

  deallog << "OK" << std::endl << std::endl;

  // make sure that is_serial_vector< dealii::PETScWrappers::MPI::BlockVector >
  // is working
  Assert(is_serial_vector<dealii::PETScWrappers::MPI::BlockVector>::value ==
           false,
         ExcInternalError());

  deallog << is_serial_vector<dealii::PETScWrappers::MPI::BlockVector>::value
          << std::endl;

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  try
    {
      test();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
