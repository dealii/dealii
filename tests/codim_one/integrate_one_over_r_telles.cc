// ---------------------------------------------------------------------
// Copyright (C) 2005 - 2018 by the deal.II authors
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

using namespace std;

int
main()
{
  initlog();
  deallog << std::fixed;

  deallog << endl
          << "Calculation of the integral of f(x,y)*1/R on [0,1]x[0,1]" << endl
          << "for f(x,y) = x^i y^j, with i,j ranging from 0 to 5, and R being"
          << endl
          << "the distance from (x,y) to four vertices of the square." << endl
          << endl;


  std::vector<Point<2>> vertices = FE_Q<2>(1).get_unit_support_points();

  for (unsigned int m = 1; m < 7; ++m)
    {
      deallog << " =========Quadrature Order: " << m
              << " =============================== " << endl;
      deallog
        << " ============================================================ "
        << endl;
      unsigned int index = 0;
      {
        deallog << " ===============Vertex: " << vertices[index]
                << " ============================= " << endl;
        QTelles<2>        quad(m, vertices[index]);
        QGaussOneOverR<2> quad2(m, vertices[index]);


        for (unsigned int i = 0; i < 6; ++i)
          {
            for (unsigned int j = 0; j < 6; ++j)
              {
                double exact_integral  = exact_integral_one_over_r(index, i, j);
                double approx_integral = 0;
                double approx_integral_2 = 0;

                for (unsigned int q = 0; q < quad.size(); ++q)
                  {
                    double x = quad.point(q)[0];
                    double y = quad.point(q)[1];
                    double R = sqrt(x * x + y * y);
                    approx_integral += (pow(x, (double)i) * pow(y, (double)j) /
                                        R * quad.weight(q));
                  }

                for (unsigned int q = 0; q < quad2.size(); ++q)
                  {
                    double x = quad2.point(q)[0];
                    double y = quad2.point(q)[1];
                    double R = sqrt(x * x + y * y);
                    approx_integral_2 +=
                      (pow(x, (double)i) * pow(y, (double)j) * quad2.weight(q));
                  }

                deallog << "f(x,y) = x^" << i << " y^" << j
                        << ", Errors = " << approx_integral - exact_integral
                        << ", " << approx_integral_2 - exact_integral
                        << std::endl;
              }
          }
      }
    }
}
