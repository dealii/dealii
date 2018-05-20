// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "interpolate_common.h"
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

//
// Check convert_generalized_support_opint_values_to_dof_values for systems
// of non-Lagrangian elements.
//

template <int dim, typename T>
void
check(T function, const unsigned int degree)
{
  FESystem<dim> fe(FE_RaviartThomas<dim>(degree),
                   2,
                   FESystem<dim>(FE_RaviartThomas<dim>(degree), 2),
                   1);
  deallog << fe.get_name() << std::endl;

  std::vector<double> dofs(fe.dofs_per_cell);

  std::vector<Vector<double>> values(fe.get_generalized_support_points().size(),
                                     Vector<double>(4 * dim));
  function.vector_value_list(fe.get_generalized_support_points(), values);

  fe.convert_generalized_support_point_values_to_dof_values(values, dofs);
  deallog << " vector " << vector_difference(fe, dofs, function, 0)
          << std::endl;
}

int
main()
{
  initlog();

  check<2>(Q1WedgeFunction<2, 1, 4 * 2>(), 1);
  check<2>(Q1WedgeFunction<2, 2, 4 * 2>(), 2);
  check<2>(Q1WedgeFunction<2, 3, 4 * 2>(), 3);
  check<3>(Q1WedgeFunction<3, 1, 4 * 3>(), 1);
  check<3>(Q1WedgeFunction<3, 2, 4 * 3>(), 2);
  check<3>(Q1WedgeFunction<3, 3, 4 * 3>(), 3);
}
