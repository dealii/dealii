// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
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

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include "interpolate_common.h"


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
