// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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


void
test()
{
  deallog << Utilities::string_to_double(" 413 ") << std::endl;

  std::vector<std::string> v;
  v.push_back("1.5");
  v.push_back(" -12.5");
  v.push_back("+125.5 ");
  AssertThrow(Utilities::string_to_double(v).size() == 3, ExcInternalError());
  deallog << Utilities::string_to_double(v)[0] << std::endl;
  deallog << Utilities::string_to_double(v)[1] << std::endl;
  deallog << Utilities::string_to_double(v)[2] << std::endl;
}



int
main()
{
  initlog();

  test();
}
