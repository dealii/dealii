// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nedelec.h>

#include "interpolate_common.h"

// FE_Nedelec<dim>::interpolate(...)

template <int dim>
void
check1(const Function<dim> &f, const unsigned int degree)
{
  FE_Nedelec<dim> fe(degree);
  deallog << fe.get_name() << ' ';
  deallog << fe.get_generalized_support_points().size() << ' ';

  std::vector<double> dofs(fe.dofs_per_cell);

  std::vector<Vector<double>> values(fe.get_generalized_support_points().size(),
                                     Vector<double>(dim));
  f.vector_value_list(fe.get_generalized_support_points(), values);
  fe.convert_generalized_support_point_values_to_dof_values(values, dofs);

  deallog << " vector " << vector_difference(fe, dofs, f, 0) << std::endl;
}

int
main()
{
  initlog();

  {
    Q1WedgeFunction<2, 0, 2> w20;
    check1(w20, 0);
    check1(w20, 1);
    check1(w20, 2);
    Q1WedgeFunction<2, 1, 2> w21;
    check1(w21, 2);
    Q1WedgeFunction<2, 2, 2> w22;
    check1(w22, 2);
    Q1WedgeFunction<2, 3, 2> w23;
    check1(w23, 3);

    PolynomialFunction<2, 1, 2> p21;
    check1(p21, 1);
    PolynomialFunction<2, 2, 2> p22;
    check1(p22, 2);
    PolynomialFunction<2, 3, 2> p23;
    check1(p23, 3);
  }

  {
    Q1WedgeFunction<3, 1, 3> w21;
    check1(w21, 1);
    check1(w21, 2);



    // FIXME - higher order interpolation is currently broken
    // PolynomialFunction<3, 1, 3> p21;
    // check1(p21, 1);
    // Q1WedgeFunction<3, 2, 3> w22;
    // check1(w22, 2);
    // Q1WedgeFunction<3, 3, 3> w23;
    // check1(w23, 3);
    // PolynomialFunction<3, 2, 3> p22;
    // check1(p22, 2);
    // PolynomialFunction<3, 3, 3> p23;
    // check1(p23, 3);
  }
}
