// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2022 by the deal.II authors
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
