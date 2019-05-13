// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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
