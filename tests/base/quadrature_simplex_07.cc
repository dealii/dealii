// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// integrates the function *f(x,y)/R, where f(x,y) is a power of x and
// y on the set [0,1]x[0,1]. dim = 2 only.

#include <deal.II/base/utilities.h>

#include "../tests.h"

// all include files needed for the program
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include "simplex.h"


int
main()
{
  initlog();

  deallog << std::endl
          << "Calculation of the integral of f(x,y)*1/R on [0,1]x[0,1]"
          << std::endl
          << "for f(x,y) = x^i y^j, with i,j ranging from 0 to 5, and R being"
          << std::endl
          << "the distance from (x,y) to [0.5,0.5]." << std::endl
          << std::endl;

  double eps = 1e-10;

  //           m  i  j
  double error[5][6][6] = {{{0}}};

  for (unsigned int m = 0; m < 5; ++m)
    {
      auto split_point = Point<2>(.5, .5);

      QSplit<2> quad(QTrianglePolar(m + 1), split_point);

      for (unsigned int i = 0; i < 6; ++i)
        for (unsigned int j = 0; j < 6; ++j)
          {
            double exact_integral  = exact_integral_one_over_r_middle(i, j);
            double approx_integral = 0;

            for (unsigned int q = 0; q < quad.size(); ++q)
              {
                double x = quad.point(q)[0];
                double y = quad.point(q)[1];
                approx_integral +=
                  (pow(x, (double)i) * pow(y, (double)j) * quad.weight(q) /
                   (quad.point(q) - split_point).norm());
              }
            error[m][i][j] = approx_integral - exact_integral;
          }
    }

  for (unsigned int i = 0; i < 6; ++i)
    for (unsigned int j = 0; j < 6; ++j)
      {
        deallog << "======= f(x,y) = x^" << i << " y^" << j << std::endl;

        for (unsigned int m = 0; m < 5; ++m)
          deallog << "Order[" << m + 1 << "], error = " << error[m][i][j]
                  << std::endl;
      }
}
