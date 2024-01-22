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
