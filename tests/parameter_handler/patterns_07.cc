// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>

#include <memory>

#include "../tests.h"

using namespace Patterns::Tools;

int
main()
{
  initlog();

  std::map<unsigned int, double> a;
  a[3] = 1.0;
  a[2] = 3.0;


  ParameterHandler prm;
  prm.add_parameter("A map", a);

  prm.log_parameters(deallog);

  prm.set("A map", "1:2.0, 3:4.0");

  deallog << "After ParameterHandler::set =========================="
          << std::endl
          << std::endl;
  prm.log_parameters(deallog);

  deallog << "Actual variables            =========================="
          << std::endl
          << std::endl;

  for (auto i : a)
    deallog << i.first << ':' << i.second << std::endl;

  return 0;
}
