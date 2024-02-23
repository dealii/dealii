// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that we set parameters for multiple instances of the same class
// correctly

#include <deal.II/base/parameter_acceptor.h>

#include "../tests.h"

struct Foo : public dealii::ParameterAcceptor
{
  Foo()
    : dealii::ParameterAcceptor("A subsection")
  {
    add_parameter("A parameter", a);
  }
  int a = 1;
};

int
main()
{
  initlog();

  Foo foo_1;
  Foo foo_2;
  Foo foo_3;

  deallog << foo_1.a << " - " << foo_2.a << " - " << foo_3.a << std::endl;

  std::stringstream parameters;
  parameters << "subsection A subsection\n"
             << "set A parameter = 2\n"
             << "end" << std::endl;
  ParameterAcceptor::initialize(parameters);

  deallog << foo_1.a << " - " << foo_2.a << " - " << foo_3.a << std::endl;
}
