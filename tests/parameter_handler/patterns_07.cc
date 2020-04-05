// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


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
    deallog << i.first << ":" << i.second << std::endl;

  return 0;
}
