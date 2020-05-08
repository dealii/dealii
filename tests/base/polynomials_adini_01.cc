// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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
      deallog << ")" << std::endl;

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
