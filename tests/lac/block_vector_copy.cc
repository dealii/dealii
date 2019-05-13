// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2018 by the deal.II authors
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



#include <deal.II/lac/block_vector.h>

#include <iostream>
#include <list>
#include <vector>

#include "../tests.h"

void
test()
{
  std::vector<double> v(9);
  for (unsigned int i = 0; i < v.size(); ++i)
    v[i] = double(i + 1);

  std::vector<types::global_dof_index> partition(3);
  for (unsigned int i = 0; i < partition.size(); ++i)
    partition[i] = 3;

  dealii::BlockVector<double> b(partition);
  AssertThrow(b.n_blocks() == partition.size(), ExcInternalError());

  unsigned int size = 0;
  for (unsigned int i = 0; i < b.n_blocks(); ++i)
    {
      AssertThrow(b.block(i).size() == partition[i], ExcInternalError());
      size += b.block(i).size();
    }
  AssertThrow(b.size() == size, ExcInternalError());

  for (unsigned int i = 0; i < b.size(); ++i)
    {
      b(i) = v[i];
      AssertThrow(b(i) == v[i], ExcInternalError());
    }

  dealii::BlockVector<double> c;
  c = b;
  AssertThrow(c == b, ExcInternalError());
  AssertThrow(c.n_blocks() == b.n_blocks(), ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  // do the same weird stuff as in
  // tests/base/reference.cc
#if __GNUC__ != 2
  std::basic_streambuf<char> *old_cerr_buf = std::cerr.rdbuf();
#else
  streambuf *old_cerr_buf = std::cerr.rdbuf();
#endif
  std::cerr.rdbuf(logfile.rdbuf());

  try
    {
      test();
    }
  catch (std::exception &e)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << e.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      // abort
      return 0;
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
      // abort
      return 0;
    };

  std::cerr.rdbuf(old_cerr_buf);
}
