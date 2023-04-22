//-----------------------------------------------------------
//
//    Copyright (C) 2023 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------


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
