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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Tests Trilinos sparsity iterators, specifically copy construction
// and assignment operators

#include "../tests.h"
#include <deal.II/lac/trilinos_sparsity_pattern.h>

void
test()
{
  TrilinosWrappers::SparsityPattern sp;

  sp.reinit(5, 7, 3);
  for(unsigned int i = 0; i < 5; ++i)
    for(unsigned int j = 0; j < 7; ++j)
      if((i + 2 * j + 1) % 3 == 0)
        {
          deallog << "Creating sparsity pattern entry " << i << ' ' << j
                  << std::endl;
          sp.add(i, j);
        }
  sp.compress();

  for(TrilinosWrappers::SparsityPattern::const_iterator p = sp.begin();
      p != sp.end();
      ++p)
    {
      deallog << p->row() << ' ' << p->column() << std::endl;

      // check copy construction
      TrilinosWrappers::SparsityPattern::const_iterator q(p);
      Assert(p == q, ExcInternalError());

      // also check copy operation
      q = p;
      Assert(p == q, ExcInternalError());
    }
}

int
main(int argc, char** argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      test();
    }
  catch(std::exception& exc)
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
  catch(...)
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
