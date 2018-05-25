// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// test add_parameters with tuples.

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/std_cxx14/memory.h>

#include <memory>

#include "../tests.h"

using namespace Patterns;
using namespace Patterns::Tools;

int
main()
{
  initlog();

  typedef std::tuple<std::string, Point<3>, unsigned int> T;

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
