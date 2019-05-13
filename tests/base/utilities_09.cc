// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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


// test Utilities::split_string_list with an empty string

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test()
{
  // test an empty string -- should yield a list of zero elements with
  // any delimiter
  {
    const char *p = "";
    deallog << Utilities::split_string_list(p).size() << std::endl;
    deallog << Utilities::split_string_list(p, ' ').size() << std::endl;
  }

  // also test a string that consists only of whitespace. this should
  // yield a list of zero elements even if (maybe not very usefully)
  // the delimiter is chosen as a whitespace itself
  {
    const char *p = "  ";
    deallog << Utilities::split_string_list(p).size() << std::endl;
    deallog << Utilities::split_string_list(p, ' ').size() << std::endl;
  }
}



int
main()
{
  initlog();

  test();
}
