// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that the FEValuesExtractors::get_name() work as expected

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include "../tests.h"

int
main()
{
  initlog();

  FEValuesExtractors::Scalar             scalar(0);
  FEValuesExtractors::Vector             vector(1);
  FEValuesExtractors::Tensor<2>          tensor(3);
  FEValuesExtractors::SymmetricTensor<2> symmetric_tensor(4);


  deallog << scalar.get_name() << std::endl
          << vector.get_name() << std::endl
          << tensor.get_name() << std::endl
          << symmetric_tensor.get_name() << std::endl;
}
