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

// test DynamicSparsityPattern::iterator with sparsity patterns that
// have an associated IndexSet

#include "../tests.h"
#include <deal.II/lac/dynamic_sparsity_pattern.h>

void
test()
{
  IndexSet rows(5);
  rows.add_index(1);
  rows.add_index(2);
  rows.add_index(4);
  DynamicSparsityPattern sp(5, 5, rows);
  sp.add(1, 1);
  sp.add(2, 2);
  sp.add(4, 4);
  sp.compress();

  DynamicSparsityPattern::const_iterator i = sp.begin();
  for(; i != sp.end(); ++i)
    deallog << i->row() << ' ' << i->column() << std::endl;

  deallog << "OK" << std::endl;

  i = sp.begin(1);
  deallog << i->row() << ' ' << i->column() << std::endl;
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
