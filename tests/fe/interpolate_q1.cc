// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include "interpolate_common.h"


// FE_Q<dim>::interpolate(...)
// FE_DGQ<dim>::interpolate(...)

template <int dim>
void
check(const Function<dim> &f, const unsigned int degree)
{
  FE_Q<dim> fe(degree);
  deallog << fe.get_name() << ' ';

  std::vector<double> dofs(fe.dofs_per_cell);

  std::vector<Vector<double>> values(fe.get_unit_support_points().size(),
                                     Vector<double>(1));
  f.vector_value_list(fe.get_unit_support_points(), values);
  fe.convert_generalized_support_point_values_to_dof_values(values, dofs);
  deallog << " vector " << difference(fe, dofs, f) << std::endl;
}

template <int dim>
void
check_dg(const Function<dim> &f, const unsigned int degree)
{
  FE_DGQ<dim> fe(degree);
  deallog << fe.get_name() << ' ';

  std::vector<double> dofs(fe.dofs_per_cell);

  std::vector<Vector<double>> values(fe.get_unit_support_points().size(),
                                     Vector<double>(1));
  f.vector_value_list(fe.get_unit_support_points(), values);
  fe.convert_generalized_support_point_values_to_dof_values(values, dofs);
  deallog << " vector " << difference(fe, dofs, f) << std::endl;
}

template <int dim>
void
check_dg_lobatto(const Function<dim> &f, const unsigned int degree)
{
  QGaussLobatto<1>          fe_quadrature(degree);
  FE_DGQArbitraryNodes<dim> fe(fe_quadrature);
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
  check_dg(w1, 1);
  check_dg(w1, 2);
  check_dg(w1, 3);
  Q1WedgeFunction<2, 1> w2;
  check(w2, 1);
  check(w2, 2);
  check(w2, 3);
  check_dg(w2, 2);
  check_dg(w2, 3);
  Q1WedgeFunction<2, 2> w22;
  check(w22, 2);
  check(w22, 3);
  check_dg(w22, 2);
  check_dg(w22, 3);
  check_dg_lobatto(w22, 4);
  Q1WedgeFunction<2, 3> w23;
  check(w23, 3);
  Q1WedgeFunction<3, 1> w3;
  check_dg(w3, 1);
  check(w3, 1);
  check(w3, 2);
  check(w3, 3);
}
