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

// test add_parameters with std::map<unsigned int,ParsedFunction>

#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/component_mask.h>

#include <memory>

#include "../tests.h"

using namespace Patterns;
using namespace Patterns::Tools;

int
main()
{
  initlog();

  typedef std::map<types::boundary_id, std::unique_ptr<FunctionParser<3>>> T;

  T a;
  a = Convert<T>::to_value("0:x,y,z*t");

  ParameterHandler prm;
  prm.add_parameter("Boundary conditions", a);

  prm.log_parameters(deallog);

  prm.set("Boundary conditions", "0:x*x*x,y*y*y,z*z*z*t;1:0,x*y,t");

  deallog << "After ParameterHandler::set =========================="
          << std::endl
          << std::endl;
  prm.log_parameters(deallog);

  deallog << "Actual variables            =========================="
          << std::endl
          << std::endl;

  deallog << Convert<T>::to_string(a) << std::endl;
}
