// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test to check that the functions in the Physics namespace compile
// with the various auto-differentiable number types:
// Transformations
//
// AD number type: Sacado DFad

#include "../tests.h"

#include "../ad_common_tests/physics_functions_03.h"

int
main()
{
  initlog();

  deallog.push("Double");
  {
    test_physics<2, double, AD::NumberTypes::sacado_dfad>();
    test_physics<3, double, AD::NumberTypes::sacado_dfad>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog.push("Float");
  {
    test_physics<2, float, AD::NumberTypes::sacado_dfad>();
    test_physics<3, float, AD::NumberTypes::sacado_dfad>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
