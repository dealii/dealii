// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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


// Evaluate FE_SimplexP_Bubbles.


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"

using namespace dealii;

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
      FE_SimplexP_Bubbles<dim, spacedim> fe(degree);
      deallog << "element tensor degree = " << fe.tensor_degree() << std::endl;
      Quadrature<dim> quad(fe.get_unit_support_points());
      test(fe, quad);
    }
}



int
main()
{
  initlog();

  test_unit_support_points<1>();
  test_unit_support_points<1, 2>();
  test_unit_support_points<1, 3>();

  test_unit_support_points<2>();
  test_unit_support_points<2, 3>();

  test_unit_support_points<3>();
}
