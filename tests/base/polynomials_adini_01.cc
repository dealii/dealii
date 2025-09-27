// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// evaluate polynomials_adini on the reference cell using a 4x4 grid

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomials_adini.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <vector>

#include "../tests.h"

template <int dim>
void
check_poly_q(const PolynomialsAdini<dim> &poly)
{
  std::vector<Point<dim>>     points;
  std::vector<double>         values;
  std::vector<Tensor<1, dim>> grads;
  std::vector<Tensor<2, dim>> grads2;
  std::vector<Tensor<3, dim>> thirds;
  std::vector<Tensor<4, dim>> fourths;

  // set up evaluation points - 4x4 points for cubic
  for (unsigned int i = 0; i < 4; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j)
        {
          Point<dim> p;
          p[0] = 1. / 3. * j;
          p[1] = 1. / 3. * i;
          points.push_back(p);
        }
    }

  // loop over evaluation points
  for (unsigned int i = 0; i < points.size(); ++i)
    {
      values.clear();
      grads.clear();
      grads2.clear();
      thirds.clear();
      fourths.clear();
      values.resize(12);
      grads.resize(12);

      deallog << "Adini<" << dim << "> point " << i << " (" << points[i][0];
      for (unsigned int d = 1; d < dim; ++d)
        deallog << ", " << points[i][d];
      deallog << ')' << std::endl;

      poly.evaluate(points[i], values, grads, grads2, thirds, fourths);

      // loop through shape fxns
      for (unsigned int j = 0; j < 12; ++j)
        {
          deallog << "Adini<" << dim << "> shape fxn " << j << ": ";
          deallog << '\t' << values[j];
          deallog << std::endl;

          deallog << "Adini<" << dim << "> shape grad " << j << ": ";
          for (unsigned int d = 0; d < dim; ++d)
            deallog << '\t' << grads[j][d];
          deallog << std::endl;
        }
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  PolynomialsAdini<2> p_2d;
  check_poly_q(p_2d);
}
