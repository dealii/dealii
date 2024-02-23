// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
