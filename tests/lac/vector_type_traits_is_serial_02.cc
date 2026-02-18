// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// check is_serial_vector type trait

#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"

void
test()
{
  // make sure that is_serial_vector< dealii::TrilinosWrappers::MPI::Vector > is
  // working
  Assert(is_serial_vector<dealii::TrilinosWrappers::MPI::Vector>::value ==
           false,
         ExcInternalError());

  deallog << is_serial_vector<dealii::TrilinosWrappers::MPI::Vector>::value
          << std::endl;

  deallog << "OK" << std::endl << std::endl;


  // make sure that dealii::TrilinosWrappers::MPI::BlockVector > is working
  Assert(is_serial_vector<dealii::TrilinosWrappers::MPI::BlockVector>::value ==
           false,
         ExcInternalError());

  deallog << is_serial_vector<dealii::TrilinosWrappers::MPI::BlockVector>::value
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
