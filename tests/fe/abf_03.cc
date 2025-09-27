// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check that convert_generalized_support_point_values_to_dof_values
// gives the correct when asking for we the nodal values for a
// function in the FE_ABF<dim(degree) ansatz space.

#include <deal.II/fe/fe_abf.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <unsigned int dim>
void
test(const unsigned int degree)
{
  FE_ABF<dim> fe(degree);

  std::vector<double> dof_values(fe.dofs_per_cell);
  for (unsigned int i = 0; i < dof_values.size(); ++i)
    dof_values[i] = 1. + 2. * random_value<double>();

  const std::vector<Point<dim>> &generalized_support_points =
    fe.get_generalized_support_points();
  std::vector<Vector<double>> real_values(generalized_support_points.size(),
                                          Vector<double>(dim));

  for (unsigned int i = 0; i < generalized_support_points.size(); ++i)
    for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
      for (unsigned int c = 0; c < dim; ++c)
        real_values[i][c] +=
          dof_values[j] *
          fe.shape_value_component(j, generalized_support_points[i], c);

  std::vector<double> compare_values(fe.dofs_per_cell);
  fe.convert_generalized_support_point_values_to_dof_values(real_values,
                                                            compare_values);

  for (unsigned int i = 0; i < dof_values.size(); ++i)
    if (std::abs(dof_values[i] - compare_values[i]) >
        std::abs(dof_values[i]) * 1.e-6)
      {
        deallog << i << ": " << dof_values[i] << " vs. " << compare_values[i]
                << std::endl;
        AssertThrow(false, ExcInternalError());
      }
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  deallog.push("2d");
  deallog.push("0");
  test<2>(0);
  deallog.pop();
  deallog.push("1");
  test<2>(1);
  deallog.pop();
  deallog.push("2");
  test<2>(2);
  deallog.pop();
  deallog.pop();

  deallog.push("3d");
  deallog.push("0");
  test<3>(0);
  deallog.pop();
  // deallog.push("1");
  // test<3>(1);
  // deallog.pop();
  deallog.pop();
}
