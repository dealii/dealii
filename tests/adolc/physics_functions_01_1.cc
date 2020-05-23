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


// Test to check that the functions in the Physics namespace compile
// with the various auto-differentiable number types:
// Kinematics
//
// AD number type: ADOL-C taped

#include "../tests.h"

#include "../ad_common_tests/physics_functions_01.h"

int
main()
{
  initlog();

  deallog.push("Double");
  {
    test_physics<2, double, AD::NumberTypes::adolc_taped>();
    test_physics<3, double, AD::NumberTypes::adolc_taped>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
