// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
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

// this is a variant of parameter_handler_15 but in fact reads data
// from a file

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
test()
{
  ParameterHandler foo;
  foo.enter_subsection("bar");
  foo.declare_entry("val", "1.0", dealii::Patterns::Double(), "");
  foo.leave_subsection();
  foo.declare_entry("val2", "2.0", dealii::Patterns::Double(), "");

  foo.parse_input(SOURCE_DIR "/parameter_handler_16_in.prm");



  foo.enter_subsection("bar");
  deallog << foo.get("val") << std::endl;
  foo.leave_subsection();
  deallog << foo.get("val2") << std::endl;
}

int
main()
{
  initlog();

  test();

  return 0;
}
