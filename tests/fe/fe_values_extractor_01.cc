// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that the FEValuesExtractors are copyable

#include <deal.II/fe/fe_values.h>

#include "../tests.h"



int
main()
{
  initlog();
  deallog << std::setprecision(2);

  {
    std::vector<FEValuesExtractors::Scalar> x;
    x.push_back(FEValuesExtractors::Scalar(42));
  }

  {
    std::vector<FEValuesExtractors::Vector> x;
    x.push_back(FEValuesExtractors::Vector(42));
  }

  deallog << "OK" << std::endl;
}
