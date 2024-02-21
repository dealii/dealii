// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/block_vector.h>

#include "../tests.h"

#define PRINTME(name, var)                              \
  deallog << "Block vector: " name << ':' << std::endl; \
  for (unsigned int i = 0; i < var.n_blocks(); ++i)     \
    deallog << "[block " << i << " ]  " << var.block(i) << std::endl;

int
main()
{
  initlog();

  {
    BlockVector<double> temp(5, 2);
    BlockVector<double> w(std::move(temp));
    PRINTME("move constructor", w);
  }

  deallog << std::endl;

  {
    BlockVector<double> u(5, 2);
    for (unsigned int i = 0; i < 5; ++i)
      for (unsigned int j = 0; j < 2; ++j)
        u.block(i)[j] = (double)(10 * i + j);

    PRINTME("BlockVector", u);

    BlockVector<double> v;
    v = u;
    PRINTME("copy assignment", v);
    PRINTME("old object", u);

    v = 0.;
    v = std::move(u);
    PRINTME("move assignment", v);
    deallog << "old object size: " << u.n_blocks() << std::endl;

    // and swap again with different sizes
    BlockVector<double> w;
    w = std::move(v);
    PRINTME("move assignment", w);
    deallog << "old object size: " << v.n_blocks() << std::endl;
  }
}
