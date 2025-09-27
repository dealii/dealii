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
