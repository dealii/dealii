// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

// test for class ArrayView

#include "../tests.h"

#include <deal.II/base/array_view.h>

template <typename T>
void
test(const T& t)
{
  try
    {
      make_array_view(std::begin(t), std::end(t));
    }
  catch(ExceptionBase& e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}

int
main()
{
  deal_II_exceptions::disable_abort_on_exception();
  initlog();

  auto d1  = std::deque<int>(2);
  auto d2  = std::deque<int>(4000);
  int  c[] = {1, 2, 3};
  auto a   = std::array<int, 3>{{1, 2, 3}};
  auto s   = std::string{"Hello world!"};
  auto v   = std::vector<int>{
    1,
    2,
    3,
  };

  deallog << "Testing std::deque<int>(2)" << std::endl;
  test(d1);
  deallog << "Testing std::deque<int>(4000)" << std::endl;
  test(d2);
  deallog << "Testing int c[3]{1,2,3}" << std::endl;
  test(c);
  deallog << "Testing std::array<int, 3> {{ 1,2,3}}" << std::endl;
  test(a);
  deallog << "Testing std::string {\"Hello world!\"}" << std::endl;
  test(s);
  deallog << "Testing std::vector<int>{1,2,3}" << std::endl;
  test(v);
}
