// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
  virtual ~FunctionBase()
  {}

  virtual void
  do_test() = 0;
};

class Function : public FunctionBase
{
public:
  Function()
    : size_(2)
  {
    deallog << "Construct object" << std::endl;
  }

  Function(const Function &f)
    : size_(f.size_)
    , vec(f.vec)
  {
    deallog << "Copy construct object" << std::endl;
  }

  virtual ~Function() override
  {
    deallog << "Destruct with size " << vec.size() << std::endl;
  }

  virtual void
  do_test() override
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
