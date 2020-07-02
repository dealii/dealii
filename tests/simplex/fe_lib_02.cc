// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Evaluate Simplex::FE_P  and Simplex::FE_DGP at quadrature points.


#include <deal.II/simplex/fe_lib.h>
#include <deal.II/simplex/quadrature_lib.h>

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
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); i++)
        deallog << fe.shape_value(i, point) << " ";
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); i++)
        deallog << fe.shape_grad(i, point) << " ";
      deallog << std::endl;
    }
  deallog << std::endl;

  deallog.pop();
}

int
main()
{
  initlog();

  test(Simplex::FE_P<2>(1), Simplex::PGauss<2>(3));
  test(Simplex::FE_P<2>(2), Simplex::PGauss<2>(7));
  test(Simplex::FE_P<3>(1), Simplex::PGauss<3>(4));
  test(Simplex::FE_P<3>(2), Simplex::PGauss<3>(10));
  test(Simplex::FE_DGP<2>(1), Simplex::PGauss<2>(3));
  test(Simplex::FE_DGP<2>(2), Simplex::PGauss<2>(7));
  test(Simplex::FE_DGP<3>(1), Simplex::PGauss<3>(4));
  test(Simplex::FE_DGP<3>(2), Simplex::PGauss<3>(10));
}
