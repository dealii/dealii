// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



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
