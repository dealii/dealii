// ------------------------------------------------------------------------
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include "../tests.h"

// integrates the function *f(x,y)/R, where f(x,y) is a power of x and
// y on the set [0,1]x[0,1]. dim = 2 only.

#include <deal.II/base/utilities.h>

// all include files needed for the program
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>

#include <string>

#include "../base/simplex.h"

// We test the integration of singular kernels with a singularity of kind 1/R
// We multiply this function with a polynomial up to degree 6.

int
main()
{
  initlog();
  deallog << std::fixed;

  deallog << std::endl
          << "Calculation of the integral of f(x,y)*1/R on [0,1]x[0,1]"
          << std::endl
          << "for f(x,y) = x^i y^j, with i,j ranging from 0 to 5, and R being"
          << std::endl
          << "the distance from (x,y) to [0.5,0.5]." << std::endl
          << std::endl;


  Point<2> center(0.5, 0.5);
  for (unsigned int m = 1; m < 6; ++m)
    {
      deallog << " =========Quadrature Order: " << m
              << " =============================== " << std::endl;
      deallog
        << " ============================================================ "
        << std::endl;
      deallog << " ===============Vertex: " << center
              << " ============================= " << std::endl;
      QTelles<2>        quad(4 * m, center);
      QGaussOneOverR<2> quad_2(m, center, true);


      for (unsigned int i = 0; i < 5; ++i)
        {
          for (unsigned int j = 0; j < 5; ++j)
            {
              double exact_integral    = exact_integral_one_over_r_middle(i, j);
              double approx_integral   = 0;
              double approx_integral_2 = 0;

              for (unsigned int q = 0; q < quad.size(); ++q)
                {
                  double x = quad.point(q)[0];
                  double y = quad.point(q)[1];
                  double R = sqrt((x - center[0]) * (x - center[0]) +
                                  (y - center[1]) * (y - center[1]));
                  approx_integral += (pow(x, (double)i) * pow(y, (double)j) /
                                      R * quad.weight(q));
                }

              for (unsigned int q = 0; q < quad_2.size(); ++q)
                {
                  double x = quad_2.point(q)[0];
                  double y = quad_2.point(q)[1];
                  double R = sqrt((x - center[0]) * (x - center[0]) +
                                  (y - center[1]) * (y - center[1]));
                  approx_integral_2 += (pow(x, (double)i) * pow(y, (double)j) /
                                        R * quad_2.weight(q));
                }


              deallog << "f(x,y) = x^" << i << " y^" << j
                      << ", Errors = " << approx_integral - exact_integral
                      << ", " << approx_integral_2 - exact_integral
                      << std::endl;
            }
        }
    }
}
