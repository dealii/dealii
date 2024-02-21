// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test QIterated for varying subdivisions.


#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"


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
