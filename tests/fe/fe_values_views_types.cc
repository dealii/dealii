// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test FEValuesViews::View<dim,spacedim,Extractor>

#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_values.h>

#include "../tests.h"

template <int dim, int spacedim, typename Extractor>
void
test(const Extractor &)
{
  typename FEValuesViews::View<dim, spacedim, Extractor>::
    template solution_value_type<double>
      t1{};

  typename FEValuesViews::View<dim, spacedim, Extractor>::
    template solution_gradient_type<double>
      t2{};

  deallog << "Test<" << Utilities::dim_string(dim, spacedim) << '>' << std::endl
          << Utilities::type_to_string(t1) << std::endl
          << Utilities::type_to_string(t2) << std::endl;
}

int
main()
{
  initlog();

  const FEValuesExtractors::Scalar    scalar(0);
  const FEValuesExtractors::Vector    vector(1);
  const FEValuesExtractors::Tensor<2> tensor(2);

  test<1, 1>(scalar);
  test<1, 1>(vector);
  test<1, 1>(tensor);

  test<1, 2>(scalar);
  test<1, 2>(vector);
  test<1, 2>(tensor);

  test<2, 2>(scalar);
  test<2, 2>(vector);
  test<2, 2>(tensor);

  test<2, 3>(scalar);
  test<2, 3>(vector);
  test<2, 3>(tensor);

  test<3, 3>(scalar);
  test<3, 3>(vector);
  test<3, 3>(tensor);
}
