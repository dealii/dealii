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
      t1;

  typename FEValuesViews::View<dim, spacedim, Extractor>::
    template solution_gradient_type<double>
      t2;

  deallog << "Test<" << Utilities::dim_string(dim, spacedim) << ">" << std::endl
          << Utilities::type_to_string(t1) << std::endl
          << Utilities::type_to_string(t2) << std::endl;
}

int
main()
{
  initlog();

  FEValuesExtractors::Scalar    scalar(0);
  FEValuesExtractors::Vector    vector(1);
  FEValuesExtractors::Tensor<2> tensor(2);

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
