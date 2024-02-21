// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that ParameterHandler::add_parameter() does not modify the
// default value

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"


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
