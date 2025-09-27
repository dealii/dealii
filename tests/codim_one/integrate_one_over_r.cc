// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2023 by the deal.II authors
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

#include <math.h>

#include <string>

#include "../base/simplex.h"

std::ofstream logfile("output");

int
main()
{
  deallog.attach(logfile);
  deallog << std::fixed;

  deallog << std::endl
          << "Calculation of the integral of f(x,y)*1/R on [0,1]x[0,1]"
          << std::endl
          << "for f(x,y) = x^i y^j, with i,j ranging from 0 to 5, and R being"
          << std::endl
          << "the distance from (x,y) to four vertices of the square."
          << std::endl
          << std::endl;

  double eps = 1e-10;

  for (unsigned int m = 1; m < 7; ++m)
    {
      deallog << " =========Quadrature Order: " << m
              << " =============================== " << std::endl;
      deallog
        << " ============================================================ "
        << std::endl;
      for (unsigned int index = 0; index < 4; ++index)
        {
          deallog << " ===============Vertex Index: " << index
                  << " ============================= " << std::endl;
          QGaussOneOverR<2> quad(m, index);
          QGaussOneOverR<2> quad_de(m, index, true);
          for (unsigned int i = 0; i < 6; ++i)
            {
              for (unsigned int j = 0; j < 6; ++j)
                {
                  double exact_integral =
                    exact_integral_one_over_r(index, i, j);
                  double approx_integral   = 0;
                  double approx_integral_2 = 0;

                  for (unsigned int q = 0; q < quad.size(); ++q)
                    {
                      double x = quad.point(q)[0];
                      double y = quad.point(q)[1];
                      approx_integral += (pow(x, (double)i) *
                                          pow(y, (double)j) * quad.weight(q));
                      approx_integral_2 +=
                        (pow(x, (double)i) * pow(y, (double)j) *
                         quad_de.weight(q) /
                         (quad_de.point(q) -
                          GeometryInfo<2>::unit_cell_vertex(index))
                           .norm());
                    }

                  deallog << "f(x,y) = x^" << i << " y^" << j
                          << ", Error = " << approx_integral - exact_integral;
                  if (std::abs(approx_integral - approx_integral_2) < eps)
                    deallog << std::endl;
                  else
                    deallog
                      << ", design: " << approx_integral_2 - exact_integral
                      << std::endl;
                }
            }
        }
    }
}
