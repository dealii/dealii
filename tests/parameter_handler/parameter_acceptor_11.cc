// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


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
