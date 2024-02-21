// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
  check("default,alsodefault",
        "default|nodefault|alsodefault",
        "set v=nodefault");

  check("  input 2  ,  have spaces  ", "have spaces|input 2", "");

  // check correct handling of space in input, default, and values:
  check("input 2",
        "have spaces|input 2",
        "set v=   input 2  ,   have spaces  ");

  check("",
        "double  spaces|input 2",
        "set v = double  spaces  ,  double  spaces");
}


int
main()
{
  initlog();

  test();

  return 0;
}
