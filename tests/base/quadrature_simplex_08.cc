// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// integrates the function *f(x,y)/R, where f(x,y) is a power of x and
// y on the set [0,1]x[0,1]. dim = 2 only.
// Compare QTrianglePolar and QLachatWatson

#include "../tests.h"
#include <deal.II/base/utilities.h>

// all include files needed for the program
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/geometry_info.h>
#include "simplex.h"

#include <iomanip>

int
main()
{
  initlog();

  deallog << std::endl
          << "Calculation of the integral of f(x,y)*1/R on [0,1]x[0,1]" << std::endl
          << "for f(x,y) = x^i y^j, with i,j ranging from 0 to 5, and R being" << std::endl
          << "the distance from (x,y) to [0.5,0.5]." << std::endl
          << std::endl;

  double eps = 1e-10;

  const unsigned int max_order = 5;

  //           m  i  j  quadtype
  double error[max_order][6][6][2] = {{{{0}}}};

  for (unsigned int m=0; m<max_order; ++m)
    {
      auto split_point = Point<2>(.5, .5);

      QSplit<2> quad(QTrianglePolar(m+1), split_point);
      QSplit<2> quad_de(QDuffy(m+1, 1.0), split_point);

      for (unsigned int i=0; i<6; ++i)
        for (unsigned int j=0; j<6; ++j)
          {
            double exact_integral  = exact_integral_one_over_r_middle(i,j);
            double approx_integral = 0;
            double approx_integral_de = 0;

            for (unsigned int q=0; q< quad.size(); ++q)
              {
                double x = quad.point(q)[0];
                double y = quad.point(q)[1];
                approx_integral += ( pow(x, (double)i) *
                                     pow(y, (double)j) *
                                     quad.weight(q) /
                                     (quad.point(q)-split_point).norm());
              }

            for (unsigned int q=0; q< quad_de.size(); ++q)
              {
                double x = quad_de.point(q)[0];
                double y = quad_de.point(q)[1];
                approx_integral_de += ( pow(x, (double)i) *
                                        pow(y, (double)j) *
                                        quad_de.weight(q) /
                                        (quad_de.point(q)-split_point).norm());
              }

            error[m][i][j][0] = approx_integral - exact_integral;
            error[m][i][j][1] = approx_integral_de - exact_integral;
          }
    }

  for (unsigned int i=0; i<6; ++i)
    for (unsigned int j=0; j<6; ++j)
      {
        deallog << "======= f(x,y) = x^" << i
                << " y^" << j << std::endl;

        for (unsigned int m=0; m<max_order; ++m)
          deallog << "Order[" << m + 1 << "], QTrianglePolar error = "
                  << std::setw(15) << error[m][i][j][0]
                  << " QLachatWatson error  = "
                  << error[m][i][j][1] << std::endl;
      }
}
