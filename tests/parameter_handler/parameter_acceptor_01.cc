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
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

// Classical way of working with parameters.

#include <deal.II/base/parameter_acceptor.h>

#include "../tests.h"

template <int dim>
class Test : public ParameterAcceptor
{
public:
  virtual void
  declare_parameters(ParameterHandler &prm)
  {
    prm.declare_entry("A double", "0.0", Patterns::Double(), "Documentation");
  };

  virtual void
  parse_parameters(ParameterHandler &prm)
  {
    deallog << "Double: " << prm.get_double("A double") << std::endl;
  };
};


int
main()
{
  initlog();
  Test<2> a;
  Test<1> b;

  ParameterHandler prm;
  a.declare_all_parameters(prm);
  prm.log_parameters(deallog);
}
