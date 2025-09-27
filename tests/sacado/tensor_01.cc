// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Make sure we can implicitly convert a tensor of AD types to a
// tensor of double. This required a bit of work under the hood
// because Sacado's AD types cannot be implicitly converted to
// doubles.

#include "../tests.h"

#include "../ad_common_tests/physics_functions_01.h"

int
main()
{
  initlog();

  using ADNumberType =
    AD::NumberTraits<double, AD::NumberTypes::sacado_dfad>::ad_type;

  const Tensor<2, 2, ADNumberType> F;
  Tensor<2, 2>                     x = F;

  deallog << "OK" << std::endl;
}
