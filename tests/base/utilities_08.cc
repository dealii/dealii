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


// Utilities::int_to_string produced wrong results with numbers larger
// than 10 digits (including padding)

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test()
{
  deallog << Utilities::int_to_string(9, 10) << std::endl;
  deallog << Utilities::int_to_string(99, 10) << std::endl;
  deallog << Utilities::int_to_string(999, 10) << std::endl;
  deallog << Utilities::int_to_string(9999, 10) << std::endl;
  deallog << Utilities::int_to_string(99999, 10) << std::endl;
  deallog << Utilities::int_to_string(999999, 10) << std::endl;
  deallog << Utilities::int_to_string(9999999, 10) << std::endl;
  deallog << Utilities::int_to_string(99999999, 10) << std::endl;
  deallog << Utilities::int_to_string(999999999, 10) << std::endl;
}



int
main()
{
  initlog();

  test();
}
