// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// ParameterHandler could not deal with missing endline at end of file
// or can it?
// http://code.google.com/p/dealii/issues/detail?id=126

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check(const char *content, double &v1, double &v2)
{
  ParameterHandler foo;
  foo.enter_subsection("bar");
  foo.declare_entry("val", "1.0", dealii::Patterns::Double(), "");
  foo.leave_subsection();
  foo.declare_entry("val2", "2.0", dealii::Patterns::Double(), "");

  std::stringstream ss(content);

  foo.parse_input(ss);



  foo.enter_subsection("bar");
  deallog << foo.get("val") << std::endl;
  v1 = foo.get_double("val");
  foo.leave_subsection();
  deallog << foo.get("val2") << std::endl;
  v2 = foo.get_double("val2");
}

void
test(std::string content)
{
  double v1, v2;
  check((content + "\n").c_str(), v1, v2);
  double v3, v4;
  check(content.c_str(), v3, v4);

  Assert(v1 == v3, ExcInternalError());
  Assert(v2 == v4, ExcInternalError());
}

int
main()
{
  initlog();

  test("subsection bar\nend");
  test("");
  test("set val2=-3");
  test("subsection bar\n set val=2\nend");

  return 0;
}
