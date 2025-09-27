// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
