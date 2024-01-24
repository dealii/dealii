// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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
