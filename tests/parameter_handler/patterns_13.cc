// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test add_parameters with tuples.

#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"

using namespace Patterns;
using namespace Patterns::Tools;

int
main()
{
  initlog();

  using T = std::tuple<std::string, Point<3>, unsigned int>;

  T a;
  a = Convert<T>::to_value("Ciao : 1.0, 2.0, 3.0 : 33");

  ParameterHandler prm;
  prm.add_parameter("A tuple", a);

  prm.log_parameters(deallog);

  prm.set("A tuple", "Mondo : 2.0, 3.0, 4.0 : 34");

  deallog << "After ParameterHandler::set =========================="
          << std::endl
          << std::endl;
  prm.log_parameters(deallog);

  deallog << "Actual variables            =========================="
          << std::endl
          << std::endl;

  deallog << Convert<T>::to_string(a) << std::endl;
}
