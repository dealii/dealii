// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for C++11 move operations on AlignedVector<unsigned int>

#include <deal.II/base/aligned_vector.h>

#include "../tests.h"


void
test()
{
  using VEC = AlignedVector<unsigned int>;

  VEC a(4, 2);

  deallog << "Size: " << a.size() << std::endl;

  VEC b = std::move(a);

  deallog << "Size of new VEC: " << b.size() << std::endl;
  deallog << "Size of moved VEC: " << a.size() << std::endl;

  deallog << "Contents of new VEC: ";
  for (unsigned int i = 0; i < b.size(); ++i)
    deallog << b[i] << ' ';
  deallog << std::endl;

  a.resize(6, 42);

  deallog << "Size of moved VEC after resize: " << a.size() << std::endl;

  deallog << "Contents of new VEC: ";
  for (unsigned int i = 0; i < a.size(); ++i)
    deallog << a[i] << ' ';
  deallog << std::endl;

  a = std::move(b);

  deallog << "Size of move-assigned VEC: " << a.size() << std::endl;
  deallog << "Size of moved VEC: " << b.size() << std::endl;
}


int
main()
{
  initlog();

  test();

  return 0;
}
