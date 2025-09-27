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


// check that Table<N,ComplicatedType> works properties (this really tests the
// data type underlying the values field of Table and that
// fill()/resize_fast() works properly).

// NOTE: The number of calls to the constructor/destructor depends on the
// actual implementation of TableBase. When that is changed, this test will
// typically fail even if the implementation of TableBase is correct. When
// adjusting the output, make sure to check this test with valgrind.

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
  dealii::Table<2, Function> table;
  table.reinit(dealii::TableIndices<2>(2, 1));
  table[1][0].do_test();
  table[0][0].do_test();
  table.reinit(dealii::TableIndices<2>(1, 1), false);
  table[0][0].do_test();
  table.reinit(dealii::TableIndices<2>(1, 1), true);
  table[0][0].do_test();
  table.reinit(0, 0);
  table.reinit(dealii::TableIndices<2>(2, 2));
  table[0][1].do_test();
}
