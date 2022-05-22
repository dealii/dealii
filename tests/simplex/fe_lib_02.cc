// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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


// Evaluate FE_SimplexP  and FE_SimplexDGP at quadrature points.


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

int
main()
{
  initlog();

  test(FE_SimplexP<2>(1), QGaussSimplex<2>(2));
  test(FE_SimplexP<2>(2), QGaussSimplex<2>(3));
  test(FE_SimplexP<3>(1), QGaussSimplex<3>(2));
  test(FE_SimplexP<3>(2), QGaussSimplex<3>(3));
  test(FE_SimplexDGP<2>(1), QGaussSimplex<2>(2));
  test(FE_SimplexDGP<2>(2), QGaussSimplex<2>(3));
  test(FE_SimplexDGP<3>(1), QGaussSimplex<3>(2));
  test(FE_SimplexDGP<3>(2), QGaussSimplex<3>(3));
}
