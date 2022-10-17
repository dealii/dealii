// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// test functions in namespace Utilities

#include <deal.II/base/utilities.h>

#include "../tests.h"

template <typename Number>
void
test()
{
  Number number = 1.23456789;
  assert(Utilities::round(number, 1) == Utilities::round(Number{1.2}, 1));
  deallog << "1 digit OK" << std::endl;
  assert(Utilities::round(number, 2) == Utilities::round(Number{1.23}, 2));
  deallog << "2 digit OK" << std::endl;
  assert(Utilities::round(number, 3) == Utilities::round(Number{1.235}, 3));
  deallog << "3 digit OK" << std::endl;
  assert(Utilities::round(number, 4) == Utilities::round(Number{1.2346}, 4));
  deallog << "4 digit OK" << std::endl;
  assert(Utilities::round(number, 5) == Utilities::round(Number{1.23457}, 5));
  deallog << "5 digit OK" << std::endl;
  assert(Utilities::round(number, 6) == Utilities::round(Number{1.234568}, 6));
  deallog << "6 digit OK" << std::endl;
  assert(Utilities::round(number, 7) == Utilities::round(Number{1.2345679}, 7));
  deallog << "7 digit OK" << std::endl;
  assert(Utilities::round(number, 8) ==
         Utilities::round(Number{1.23456789}, 8));
  deallog << "8 digit OK" << std::endl;
}



int
main()
{
  initlog();

  deallog << "Double" << std::endl;
  test<double>();
  deallog << "Float" << std::endl;
  test<float>();
  return 0;
}
