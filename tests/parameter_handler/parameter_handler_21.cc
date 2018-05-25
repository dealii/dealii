// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
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



// check the Patterns::MultipleSelection

#include <deal.II/base/parameter_handler.h>

#include <sstream>

#include "../tests.h"

void
check(const char *defaults, const char *defined, const char *input)
{
  ParameterHandler prm;
  prm.declare_entry("v", defaults, Patterns::MultipleSelection(defined), "");

  std::stringstream in(input);
  prm.parse_input(in);

  deallog << "defaults='" << defaults << "' defined='" << defined << "' input='"
          << input << "' result='" << prm.get("v") << "'" << std::endl;
}

void
test()
{
  check("", "one option", "set v=");
  check("one option", "one option", "");
  check("", "one option", "set v=one option");

  check("", "bla|bla 2|1", "");
  check("", "bla|bla 2|1", "set v=bla 2");
  check("", "bla|bla 2|1", "set v=1,bla 2");
  check("", "bla|bla 2|1", "set v=bla,bla,bla");

  check("default,alsodefault", "default|nodefault|alsodefault", "");
  check(
    "default,alsodefault", "default|nodefault|alsodefault", "set v=nodefault");

  check("  input 2  ,  have spaces  ", "have spaces|input 2", "");

  // check correct handling of space in input, default, and values:
  check(
    "input 2", "have spaces|input 2", "set v=   input 2  ,   have spaces  ");

  check(
    "", "double  spaces|input 2", "set v = double  spaces  ,  double  spaces");
}


int
main()
{
  initlog();

  test();

  return 0;
}
