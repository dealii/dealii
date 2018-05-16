//-----------------------------------------------------------
//
//    Copyright (C) 2017 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal.II distribution.
//
//-----------------------------------------------------------


#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/parameter_acceptor.h>

#include <deal.II/base/point.h>

template<int dim>
class Test : public ParameterAcceptor
{
public:
  Test()
  {
    std::string def = "0.";
    for (int i=1; i<dim; ++i)
      def += ",0.";
    add_parameter("A point", p);
  };

  void
  log_info()
  {
    deallog << "My type: " << boost::core::demangle(typeid(*this).name()) << std::endl
            << "p: " << p << std::endl;
  }

private:
  Point<dim> p;
};


int
main ()
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
