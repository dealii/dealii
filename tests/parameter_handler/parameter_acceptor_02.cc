// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/utilities.h>

#include <boost/core/demangle.hpp>

#include "../tests.h"

template <int dim>
class Test : public ParameterAcceptor
{
public:
  Test()
  {
    add_parameter("A double", a);
    add_parameter("An int", b);
    add_parameter("A string", c);
    add_parameter("A bool", d);
  };

  void
  log_info()
  {
    deallog << "My type: " << boost::core::demangle(typeid(*this).name())
            << std::endl
            << "a: " << a << std::endl
            << "b: " << b << std::endl
            << "c: " << c << std::endl
            << "d: " << (d ? "true" : "false") << std::endl;
  }

private:
  double      a = 1.0;
  int         b = 2;
  std::string c = "Ciao";
  bool        d = true;
};


int
main()
{
  initlog();
  Test<1> a;
  Test<2> b;

  auto &prm = ParameterAcceptor::prm;

  ParameterAcceptor::declare_all_parameters();
  ParameterAcceptor::parse_all_parameters();
  prm.log_parameters(deallog);

  a.log_info();
  b.log_info();

  prm.parse_input_from_string(""
                              "subsection Test<1>\n"
                              "  set A double = 3.0\n"
                              "end\n"
                              "subsection Test<2>\n"
                              "  set A double = 5.0\n"
                              "end\n");

  prm.log_parameters(deallog);
  ParameterAcceptor::parse_all_parameters(prm);

  a.log_info();
  b.log_info();
}
