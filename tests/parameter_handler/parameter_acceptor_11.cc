// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that clear works as expected.

#include <deal.II/base/parameter_acceptor.h>

#include "../tests.h"

struct Foo : public dealii::ParameterAcceptor
{};


struct Bar : public dealii::ParameterAcceptor
{
  Bar()
  {
    add_parameter("A parameter", a);
  }
  int a = 1;
};

int
main()
{
  initlog();
  {
    Foo foo;
    ParameterAcceptor::clear();
    // <-- foo goes out of scope here
  }

  {
    Bar bar;
    ParameterAcceptor::prm.log_parameters(deallog);
    ParameterAcceptor::clear();
    ParameterAcceptor::prm.log_parameters(deallog);
  }
}
