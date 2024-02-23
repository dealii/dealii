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


// Evaluation of a single component (scalar) system using a helper class
//
// AD number type: ADOL-C tapeless

#include "../tests.h"

#include "../ad_common_tests/helper_scalar_single_component_01.h"

int
main()
{
  initlog();

  deallog.push("Double");
  {
    test_scalar<2, double, AD::NumberTypes::adolc_tapeless>();
    test_scalar<3, double, AD::NumberTypes::adolc_tapeless>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
