// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for AlignedVector::push_back

#include <deal.II/base/aligned_vector.h>

#include "../tests.h"

// dummy class without default constructor
class foo
{
public:
  foo(const unsigned int a)
    : vec(1, a)
  {}
  foo(const foo &bar)
    : vec(bar.vec)
  {}
  foo &
  operator=(const foo &bar)
  {
    vec = bar.vec;
    return *this;
  }

  unsigned int
  element() const
  {
    return vec[0];
  }

private:
  AlignedVector<unsigned int> vec;
};

void
test()
{
  AlignedVector<foo> vec;
  vec.push_back(foo(3));
  vec.resize(2, foo(6));
  deallog << vec[0].element() << ' ' << vec[1].element() << std::endl;
  vec[1] = foo(13);
  deallog << vec[0].element() << ' ' << vec[1].element() << std::endl;
  vec.reserve(3);
  deallog << vec[0].element() << ' ' << vec[1].element() << std::endl;
}



int
main()
{
  initlog();

  test();
}
