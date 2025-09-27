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
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include "../tests.h"

template <int dim>
class Test : public ParameterAcceptor
{
public:
  Test()
  {
    std::string def = "0.";
    for (int i = 1; i < dim; ++i)
      def += ",0.";
    add_parameter("A point", p);
  };

  void
  log_info()
  {
    deallog << "My type: " << boost::core::demangle(typeid(*this).name())
            << std::endl
            << "p: " << p << std::endl;
  }

private:
  Point<dim> p;
};


int
main()
{
  initlog();
  Test<1> a;
  Test<2> b;
  Test<3> c;

  auto &prm = ParameterAcceptor::prm;
  ParameterAcceptor::declare_all_parameters();
  prm.parse_input_from_string(""
                              "subsection Test<1>\n"
                              "  set A point = 1.0\n"
                              "end\n"
                              "subsection Test<2>\n"
                              "  set A point = 1.0, 2.0\n"
                              "end\n"
                              "subsection Test<3>\n"
                              "  set A point = 1.0, 2.0, 3.0\n"
                              "end\n");

  prm.log_parameters(deallog);
  ParameterAcceptor::parse_all_parameters(prm);

  a.log_info();
  b.log_info();
  c.log_info();
}
