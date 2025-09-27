// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Tests Trilinos sparsity iterators, specifically copy construction
// and assignment operators
//
// The accessors keep a cache that is administered via a shared
// pointer. make sure that copying these objects works as expected

#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  TrilinosWrappers::SparsityPattern sp;

  sp.reinit(5, 7, 3);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 7; ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        {
          deallog << "Creating sparsity pattern entry " << i << ' ' << j
                  << std::endl;
          sp.add(i, j);
        }
  sp.compress();

  // get an iterator to the first row, then copy it
  TrilinosWrappers::SparsityPattern::const_iterator p = sp.begin(0);
  TrilinosWrappers::SparsityPattern::const_iterator q(p);

  // now reset 'p' to the second row, and output what elements both
  // 'p' and 'q' point to
  p = sp.begin(1);

  for (unsigned int i = 0; i < 2; ++i, ++p)
    deallog << "p[" << i << "]: " << p->row() << ' ' << p->column()
            << std::endl;

  for (unsigned int i = 0; i < 2; ++i, ++q)
    deallog << "q[" << i << "]: " << q->row() << ' ' << q->column()
            << std::endl;


  // do the same with copy assignment
  p = sp.begin(0);
  q = p;
  p = sp.begin(1);

  for (unsigned int i = 0; i < 2; ++i, ++p)
    deallog << "p[" << i << "]: " << p->row() << ' ' << p->column()
            << std::endl;

  for (unsigned int i = 0; i < 2; ++i, ++q)
    deallog << "q[" << i << "]: " << q->row() << ' ' << q->column()
            << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

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
