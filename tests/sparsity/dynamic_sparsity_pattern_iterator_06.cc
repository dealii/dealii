// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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



// investigate performance for DynamicSparsityPattern::begin(r) if the
// SP is empty

#include "../tests.h"
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <fstream>
#include <iomanip>


void test (bool have_set)
{
  const int size = 100000000;
  const int my_start = size/3;

  IndexSet empty_set;
  IndexSet owned(size);
  owned.add_range(my_start, my_start+5);

  DynamicSparsityPattern sp (size, 5, have_set?owned:empty_set);

  for (unsigned int i=my_start; i<my_start+5; ++i)
    for (DynamicSparsityPattern::iterator p=sp.begin(i);
         p != sp.end(i); ++p)
      deallog << p->row() << ' ' << p->column() << std::endl;

  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  try
    {
      test (false);
      test (true);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
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
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
