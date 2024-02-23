// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Ensure that FE_SimplexP_Bubbles works with FESystem at higher-order. We
// previously hardcoded some assumptions which prevented non-hypercube
// elements with more than zero DoF per quad from working with FESystem.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test(const FiniteElement<dim, spacedim> &fe, const Quadrature<dim> &quad)
{
  deallog.push(fe.get_name());

  for (const auto &point : quad.get_points())
    {
      deallog << point << " : " << std::endl;
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        deallog << fe.shape_value(i, point) << ' ';
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        deallog << fe.shape_grad(i, point) << ' ';
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        deallog << fe.shape_grad_grad(i, point) << ' ';
      deallog << std::endl;
    }
  deallog << std::endl;

  deallog.pop();
}



template <int dim, int spacedim = dim>
void
test_unit_support_points()
{
  deallog << "Test support points for dim = " << dim
          << " and spacedim = " << spacedim << std::endl;
  for (unsigned int degree = 1; degree < 3; ++degree)
    {
      deallog << "approximation degree = " << degree << std::endl;
      FESystem<dim, spacedim> fe(FE_SimplexP_Bubbles<dim, spacedim>(degree),
                                 dim);
      deallog << "element tensor degree = " << fe.tensor_degree() << std::endl;
      Quadrature<dim> quad(
        fe.reference_cell().template get_midpoint_quadrature<dim>());
      test(fe, quad);
    }
}



int
main()
{
  initlog();

  test_unit_support_points<3>();
}
