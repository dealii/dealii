// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Evaluation of a single component (tensor) system using a helper class
//
// AD number type: Sacado DFad-DFad

#include "../tests.h"

#include "../ad_common_tests/helper_scalar_single_component_03.h"

int
main()
{
  initlog();

  deallog.push("Double");
  {
    test_tensor<2, double, AD::NumberTypes::sacado_dfad_dfad>();
    test_tensor<3, double, AD::NumberTypes::sacado_dfad_dfad>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog.push("Float");
  {
    test_tensor<2, float, AD::NumberTypes::sacado_dfad_dfad>();
    test_tensor<3, float, AD::NumberTypes::sacado_dfad_dfad>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
