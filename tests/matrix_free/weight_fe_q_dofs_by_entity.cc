// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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



// Test weight_fe_q_dofs_by_entity() and
// compute_weights_fe_q_dofs_by_entity_shifted()

#include <deal.II/matrix_free/tensor_product_kernels.h>

#include <iomanip>
#include <iostream>

#include "../tests.h"

int
main()
{
  initlog();

  {
    const unsigned int dim       = 2;
    const unsigned int fe_degree = 5;
    using Number                 = double;

    std::vector<Number> weights(Utilities::pow(3, dim));
    for (unsigned int i = 0; i < weights.size(); ++i)
      weights[i] = i;

    std::vector<Number> values(Utilities::pow(fe_degree + 1, dim), 1.0);

    internal::weight_fe_q_dofs_by_entity<dim, -1, Number>(weights.data(),
                                                          1,
                                                          fe_degree + 1,
                                                          values.data());

    for (unsigned int i_1 = 0, c = 0; i_1 < fe_degree + 1; ++i_1)
      {
        for (unsigned int i_0 = 0; i_0 < fe_degree + 1; ++i_0, ++c)
          deallog << values[c] << " ";
        deallog << std::endl;
      }
    deallog << std::endl;

    for (auto &i : weights)
      i = 0.0;

    internal::compute_weights_fe_q_dofs_by_entity<dim, -1, Number>(
      values.data(), 1, fe_degree + 1, weights.data());

    for (const auto i : weights)
      deallog << i << " ";

    deallog << std::endl;
    deallog << std::endl;
  }

  {
    const unsigned int dim       = 2;
    const unsigned int fe_degree = 4;
    using Number                 = double;

    std::vector<Number> weights(Utilities::pow(3, dim));
    for (unsigned int i = 0; i < weights.size(); ++i)
      weights[i] = i;

    std::vector<Number> values(Utilities::pow((2 * fe_degree - 1), dim), 1.0);

    internal::weight_fe_q_dofs_by_entity_shifted<dim, -1, Number>(
      weights.data(), 1, 2 * fe_degree - 1, values.data());

    for (unsigned int i_1 = 0, c = 0; i_1 < (2 * fe_degree - 1); ++i_1)
      {
        for (unsigned int i_0 = 0; i_0 < (2 * fe_degree - 1); ++i_0, ++c)
          deallog << values[c] << " ";
        deallog << std::endl;
      }
    deallog << std::endl;

    for (auto &i : weights)
      i = 0.0;

    internal::compute_weights_fe_q_dofs_by_entity_shifted<dim, -1, Number>(
      values.data(), 1, 2 * fe_degree - 1, weights.data());

    for (const auto i : weights)
      deallog << i << " ";

    deallog << std::endl;
    deallog << std::endl;
  }
}
