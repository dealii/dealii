// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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

// check FiniteSizeHistory class


#include <deal.II/numerics/history.h>

#include "../tests.h"

void
test()
{
  const unsigned int        max_size = 3;
  FiniteSizeHistory<double> storage(max_size);
  deallog << "size:     " << storage.size() << std::endl
          << "max_size: " << storage.max_size() << std::endl;

  for (unsigned int i = 0; i < 2; ++i)
    storage.add(0.1 * (i + 1));

  // 2 elements
  deallog << "initial:" << std::endl;
  for (unsigned int i = 0; i < storage.size(); ++i)
    deallog << storage[i] << std::endl;

  // 3 elements
  deallog << "append:" << std::endl;
  storage.add(0.555);
  for (unsigned int i = 0; i < storage.size(); ++i)
    deallog << storage[i] << std::endl;

  // 2 elements:
  deallog << "remove second element:" << std::endl;
  storage.remove(1);
  for (unsigned int i = 0; i < storage.size(); ++i)
    deallog << storage[i] << std::endl;

  // 2 elements:
  deallog << "change 0th:" << std::endl;
  storage[0] = 0.33;
  for (unsigned int i = 0; i < storage.size(); ++i)
    deallog << storage[i] << std::endl;

  // 3 elements:
  deallog << "append:" << std::endl;
  storage.add(0.666);
  for (unsigned int i = 0; i < storage.size(); ++i)
    deallog << storage[i] << std::endl;

  // 3 elements:
  deallog << "change 0th:" << std::endl;
  storage[0] = 0.22;
  for (unsigned int i = 0; i < storage.size(); ++i)
    deallog << storage[i] << std::endl;

  // 3 elements:
  deallog << "append:" << std::endl;
  storage.add(0.777);
  for (unsigned int i = 0; i < storage.size(); ++i)
    deallog << storage[i] << std::endl;

  // 2 elements:
  deallog << "remove last:" << std::endl;
  storage.remove(storage.size() - 1);
  for (unsigned int i = 0; i < storage.size(); ++i)
    deallog << storage[i] << std::endl;

  storage.clear();
}


int
main(int argc, char **argv)
{
  initlog();

  test();
}
