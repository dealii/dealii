// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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
  catch (std::exception &exc)
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
