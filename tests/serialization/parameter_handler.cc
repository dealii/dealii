// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check serialization for Parameter_handler

#include <deal.II/base/parameter_handler.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"


void
test()
{
  ParameterHandler prm1;
  prm1.declare_entry("int1", "1", Patterns::Integer(), "doc 1");
  prm1.declare_entry("int2", "2", Patterns::Integer(), "doc 2");
  prm1.enter_subsection("ss1");
  {
    prm1.declare_entry("double 1", "1.234", Patterns::Double(), "doc 3");

    prm1.enter_subsection("ss2");
    {
      prm1.declare_entry("double 2", "4.321", Patterns::Double(), "doc 4");
    }
    prm1.leave_subsection();
  }
  prm1.leave_subsection();

  // things with strange characters
  prm1.enter_subsection("Testing%testing");
  {
    prm1.declare_entry("string&list",
                       "< & > ; /",
                       Patterns::Anything(),
                       "docs 1");
    prm1.declare_entry("int*int", "2", Patterns::Integer());
    prm1.declare_entry("double+double",
                       "6.1415926",
                       Patterns::Double(),
                       "docs 3");
  }
  prm1.leave_subsection();

  ParameterHandler prm2;
  prm2.enter_subsection("s234s1");
  {
    prm2.declare_entry("dummy 1", "19.4", Patterns::Double(), "paper 4");

    prm2.enter_subsection("s2skjds2");
    {
      prm2.declare_entry("var 2", "6.321", Patterns::Double(), "invalid");
    }
    prm2.leave_subsection();
  }
  prm2.leave_subsection();

  // things with strange characters
  prm2.enter_subsection("Testing%testing");
  {
    prm2.declare_entry("int*int", "2", Patterns::Integer());
    prm2.declare_entry("string&list",
                       "< & > ; /",
                       Patterns::Anything(),
                       "docs 1");
  }
  prm2.leave_subsection();

  prm2.declare_entry("int1", "1.", Patterns::Double(), "doc 1");
  prm2.declare_entry("int2", "2", Patterns::Anything(), "doc 2");

  ParameterHandler prm3;

  verify(prm1, prm2);

  verify(prm1, prm3);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
