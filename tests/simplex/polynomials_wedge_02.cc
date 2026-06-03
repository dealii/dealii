// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check if PolynomialsWedge is 1 at the corresponding support point and 0
// everywhere else


#include <deal.II/base/polynomials_wedge.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int degree)
{
  FE_WedgeP<dim> fe_wedge(degree);
  const auto     support_points = fe_wedge.get_unit_support_points();
  const auto     poly = ScalarLagrangePolynomialWedge(degree, support_points);

  for (unsigned int i = 0; i < support_points.size(); ++i)
    for (unsigned int j = 0; j < fe_wedge.n_dofs_per_cell(); ++j)
      {
        const auto value = poly.compute_value(j, support_points[i]);
        if (i == j)
          {
            // has to be 1
            if (std::abs(value - 1.0) < 1e-12)
              deallog << "ok ";
            else
              deallog << "Failure for shape function " << j
                      << " on support point " << i << std::endl;
          }
        else
          {
            // has to be 0
            if (std::abs(value) < 1e-12)
              deallog << "ok ";
            else
              deallog << "Failure for shape function " << j
                      << " on support point " << i << std::endl;
          }
        deallog << std::endl;
      }
}


int
main()
{
  initlog();
  for (unsigned int i = 1; i < 3; ++i)
    {
      deallog.push("3d-" + std::to_string(i));
      test<3>(i);
      deallog.pop();
    }
}
