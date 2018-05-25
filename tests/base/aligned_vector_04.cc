// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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


// check that AlignedVector::fill() does correctly call the destructor and
// constructor on a complicated class

// NOTE: The number of calls to the constructor/destructor depends on the
// actual implementation of AlignedVector. When that is changed, this test
// will typically fail even if the implementation in AlignedVector is
// otherwise correct. When adjusting the output, make sure to check this test
// with valgrind.

#include <deal.II/base/table.h>

#include "../tests.h"

// make function virtual to ensure that the function table is correctly copied
class FunctionBase
{
public:
  ~FunctionBase()
  {}

  virtual void
  do_test() = 0;
};

class Function
{
public:
  Function() : size_(2)
  {
    deallog << "Construct object" << std::endl;
  }

  Function(const Function &f) : size_(f.size_), vec(f.vec)
  {
    deallog << "Copy construct object" << std::endl;
  }

  ~Function()
  {
    deallog << "Destruct with size " << vec.size() << std::endl;
  }

  virtual void
  do_test()
  {
    vec.resize(size_++);
    deallog << "Resize vector to " << vec.size() << std::endl;
  }

private:
  unsigned int              size_;
  std::vector<unsigned int> vec;
};

int
main()
{
  initlog();
  AlignedVector<Function> vec;
  vec.resize(2);
  vec[1].do_test();
  vec[0].do_test();
  vec[0].do_test();
  vec.fill(Function());
  vec.resize(1);
  vec[0].do_test();
  vec.resize(3);
  vec[0].do_test();
  vec.fill(Function());
  vec[0].do_test();
}
