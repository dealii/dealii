// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test various Quadrature::initialize functions


#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

template <int dim>
void
compare_quadratures(const Quadrature<dim> &reference,
                    const Quadrature<dim> &check)
{
  AssertDimension(reference.size(), check.size());
  for (unsigned int i = 0; i < reference.size(); ++i)
    {
      // Check for exact equality because data gets copied, so there should be
      // no roundoff
      if (reference.point(i).distance(check.point(i)) > 0)
        deallog << " error point " << i << reference.point(i) << " vs "
                << check.point(i);
      if (reference.weight(i) != check.weight(i))
        deallog << " error weight " << i << reference.weight(i) << " vs "
                << check.weight(i);
    }
  deallog << "OK" << std::endl;
}

template <int dim>
void
test()
{
  deallog.push(std::to_string(dim) + "d");
  Quadrature<dim> quadrature_2;
  {
    QGauss<dim> quad(2);
    quadrature_2.initialize(quad.get_points(), quad.get_weights());
    deallog << "Testing initialize(vector&, vector&): ";
    compare_quadratures(quad, quadrature_2);
  }
  {
    QGauss<dim> quad(3);
    quadrature_2.initialize(quad.get_points(), quad.get_weights());
    deallog << "Testing initialize(vector&, vector&): ";
    compare_quadratures(quad, quadrature_2);
  }
  {
    QGauss<dim> quad(2);
    quadrature_2.initialize(make_array_view(quad.get_points()),
                            make_array_view(quad.get_weights()));
    deallog << "Testing initialize(ArrayView, ArrayView): ";
    compare_quadratures(quad, quadrature_2);
  }
  {
    Quadrature<dim> quad(quadrature_2.get_points());
    quadrature_2.initialize(make_array_view(quad.get_points()));
    deallog << "Testing initialize(ArrayView): ";
    compare_quadratures(quad, quadrature_2);
  }
  deallog.pop();
}

int
main()
{
  initlog();

  test<1>();
  test<2>();
}
