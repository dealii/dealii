// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check to_string and to_value

#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"

using dealii::Patterns::Tools::to_string;
using dealii::Patterns::Tools::to_value;

int
main()
{
  initlog();

  auto a = std::make_tuple(1, std::string("ciao"));

  auto s = to_string(a);
  to_value("2 : mondo", a);

  deallog << "From: " << s << " to " << to_string(a) << std::endl;
}
