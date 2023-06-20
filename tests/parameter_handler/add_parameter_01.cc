// ---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2023 by the deal.II authors
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



// check that ParameterHandler::add_parameter() does not modify the
// default value

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

using namespace dealii;

int
main()
{
  initlog();

  double a = std::numeric_limits<double>::lowest();

  AssertThrow(a == std::numeric_limits<double>::lowest(), ExcInternalError());

  ParameterHandler prm;
  prm.add_parameter("test", a);

  AssertThrow(a == std::numeric_limits<double>::lowest(), ExcInternalError());

  deallog << "OK!" << std::endl;
}
