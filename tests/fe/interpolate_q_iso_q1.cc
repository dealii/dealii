// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2020 by the deal.II authors
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

#include <deal.II/fe/fe_q_iso_q1.h>

#include "interpolate_common.h"


// FE_Q_iso_Q1<dim>::interpolate(...)

template <int dim>
void
check(const Function<dim> &f, const unsigned int degree)
{
  FE_Q_iso_Q1<dim> fe(degree);
  deallog << fe.get_name() << ' ';

  std::vector<double> dofs(fe.dofs_per_cell);

  std::vector<Vector<double>> values(fe.get_unit_support_points().size(),
                                     Vector<double>(1));
  f.vector_value_list(fe.get_unit_support_points(), values);
  fe.convert_generalized_support_point_values_to_dof_values(values, dofs);
  deallog << " vector " << difference(fe, dofs, f) << std::endl;
}

int
main()
{
  initlog();

  Q1WedgeFunction<1, 1> w1;
  check(w1, 1);
  check(w1, 2);
  check(w1, 3);
  Q1WedgeFunction<2, 1> w2;
  check(w2, 1);
  check(w2, 2);
  check(w2, 3);
  Q1WedgeFunction<3, 1> w3;
  check(w3, 1);
  check(w3, 2);
  check(w3, 3);
}
