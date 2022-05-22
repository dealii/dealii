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


// Test QIterated for varying subdivisions.


#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test(const Quadrature<1> &q_outer, const Quadrature<1> &q_inner)
{
  QIterated<dim> quad(q_inner, q_outer.get_points());

  for (unsigned int q = 0; q < quad.size(); ++q)
    {
      deallog << quad.point(q) << ' ';
      deallog << quad.weight(q) << ' ';
      deallog << std::endl;
    }
  deallog << std::endl;
}

template <int dim>
void
test()
{
  for (unsigned int i = 2; i <= 4; ++i)
    for (unsigned int j = 1; j <= 3; ++j)
      {
        deallog.push("QGaussLobatto<1>(" + std::to_string(i) +
                     ") + QGauss<1>(" + std::to_string(j) + ")");
        test<dim>(QGaussLobatto<1>(i), QGauss<1>(j));
        deallog.pop();
      }

  for (unsigned int i = 2; i <= 4; ++i)
    for (unsigned int j = 2; j <= 3; ++j)
      {
        deallog.push("QGaussLobatto<1>(" + std::to_string(i) +
                     ") + QGaussLobatto<1>(" + std::to_string(j) + ")");
        test<dim>(QGaussLobatto<1>(i), QGaussLobatto<1>(j));
        deallog.pop();
      }
}

int
main()
{
  initlog();

  test<1>();
  test<2>();
}
