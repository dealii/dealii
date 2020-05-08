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

#include "interpolate_common.h"


// FE_Q_Bubbles<dim>::interpolate(...)

template <int dim>
void
check_q_bubbles(const Function<dim> &f, const unsigned int degree)
{
  FE_Q_Bubbles<dim> fe(degree);
  deallog << fe.get_name() << ' ';

  std::vector<double> dofs(fe.dofs_per_cell);

  std::vector<Vector<double>> values(fe.get_unit_support_points().size(),
                                     Vector<double>(1));
  f.vector_value_list(fe.get_unit_support_points(), values);
  fe.convert_generalized_support_point_values_to_dof_values(values, dofs);
  deallog << " value " << difference(fe, dofs, f) << std::endl;
}

int
main()
{
  initlog();

  Q1WedgeFunction<1, 1> w1;
  check_q_bubbles(w1, 1);
  check_q_bubbles(w1, 2);
  check_q_bubbles(w1, 3);
  Q1WedgeFunction<2, 1> w2;
  check_q_bubbles(w2, 1);
  check_q_bubbles(w2, 2);
  check_q_bubbles(w2, 3);
  Q1WedgeFunction<2, 2> w22;
  check_q_bubbles(w22, 2);
  check_q_bubbles(w22, 3);
  Q1WedgeFunction<2, 3> w23;
  check_q_bubbles(w23, 3);
  Q1WedgeFunction<3, 1> w3;
  check_q_bubbles(w3, 1);
  check_q_bubbles(w3, 2);
  check_q_bubbles(w3, 3);
}
