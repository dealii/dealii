// ---------------------------------------------------------------------
// Copyright (C) 2005 - 2015 by the deal.II authors
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

#include "../tests.h"

// integrates the function *f(x,y)/R, where f(x,y) is a power of x and
// y on the set [0,1]x[0,1]. dim = 2 only.

#include <deal.II/base/utilities.h>

// all include files needed for the program
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/fe/fe_q.h>

#include <string>

using namespace std;
using namespace dealii;

// We test the integration of singular kernels with a singularity of kind 1/R
// We multiply this function with a polynomial up to degree 6.

double
exact_integral_one_over_r_middle (
  const unsigned int i, const unsigned int j);

int
main ()
{
  initlog();
  deallog << std::fixed;

  deallog << endl
          << "Calculation of the integral of f(x,y)*1/R on [0,1]x[0,1]" << endl
          << "for f(x,y) = x^i y^j, with i,j ranging from 0 to 5, and R being" << endl
          << "the distance from (x,y) to [0.5,0.5]." << endl
          << endl;


  Point <2> center(0.5,0.5);
  for (unsigned int m=1; m<6; ++m)
    {
      deallog << " =========Quadrature Order: " << m
              << " =============================== " << endl;
      deallog << " ============================================================ " << endl;
      deallog << " ===============Vertex: " << center
              << " ============================= " << endl;
      QTelles<2> quad(4*m, center);
      QGaussOneOverR<2> quad_2(m, center, true);


      for (unsigned int i = 0; i < 5; ++i)
        {
          for (unsigned int j = 0; j < 5; ++j)
            {

              double exact_integral  = exact_integral_one_over_r_middle(i,j);
              double approx_integral = 0;
              double approx_integral_2 = 0;

              for (unsigned int q=0; q< quad.size(); ++q)
                {
                  double x = quad.point(q)[0];
                  double y = quad.point(q)[1];
                  double R = sqrt((x-center[0])*(x-center[0])+(y-center[1])*(y-center[1]));
                  approx_integral += ( pow(x, (double)i) *
                                       pow(y, (double)j) / R *
                                       quad.weight(q) );
                }

              for (unsigned int q=0; q< quad_2.size(); ++q)
                {
                  double x = quad_2.point(q)[0];
                  double y = quad_2.point(q)[1];
                  double R = sqrt((x-center[0])*(x-center[0])+(y-center[1])*(y-center[1]));
                  approx_integral_2 += ( pow(x, (double)i) *
                                         pow(y, (double)j) / R *
                                         quad_2.weight(q) );
                }


              deallog << "f(x,y) = x^" << i
                      << " y^" << j
                      << ", Errors = "
                      << approx_integral - exact_integral
                      << ", "
                      << approx_integral_2 - exact_integral
                      << std::endl;
            }
        }
    }
}

double exact_integral_one_over_r_middle(const unsigned int i,
                                        const unsigned int j)
{
  Assert(i<6, ExcNotImplemented());
  Assert(j<6, ExcNotImplemented());

// The integrals are computed using the following Mathematica snippet of
// code:
//
// x0 = 0.5
// y0 = 0.5
// Do[Do[Print["v[0][", n, "][", m, "]=",
//    NumberForm[
//     NIntegrate[
//      x^n*y^m/Sqrt[(x - x0)^2 + (y - y0)^2], {x, 0, 1}, {y, 0, 1},
//      MaxRecursion -> 10000, PrecisionGoal -> 9], 9], ";"], {n, 0,
//    4}], {m, 0, 4}]


  static double v[1][6][6] =
  {
    {
      { 0}
    }
  };
  if (v[0][0][0] == 0)
    {
      v[0][0][0] = 3.52549435;;
      v[0][1][0] = 1.76274717;
      v[0][2][0]=1.07267252;
      v[0][3][0]=0.727635187;
      v[0][4][0]=0.53316959;
      v[0][0][1]=1.76274717;
      v[0][1][1]=0.881373587;
      v[0][2][1]=0.536336258;
      v[0][3][1]=0.363817594;
      v[0][4][1]=0.266584795;
      v[0][0][2]=1.07267252;
      v[0][1][2]=0.536336258;
      v[0][2][2]=0.329313861;
      v[0][3][2]=0.225802662;
      v[0][4][2]=0.167105787;
      v[0][0][3]=0.727635187;
      v[0][1][3]=0.363817594;
      v[0][2][3]=0.225802662;
      v[0][3][3]=0.156795196;
      v[0][4][3]=0.117366283;
      v[0][0][4]=0.53316959;
      v[0][1][4]=0.266584795;
      v[0][2][4]=0.167105787;
      v[0][3][4]=0.117366283;
      v[0][4][4]=0.0887410133;
    }
  return v[0][i][j];
}
