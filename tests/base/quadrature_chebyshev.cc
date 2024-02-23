// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check accuracy of the Chebyshev quadrature formulas by using them to
// integrate polynomials of increasing degree, and finding the degree
// until which they integrate exactly


#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"



template <typename quadrature_type, unsigned short startn>
void
check_quadrature(double *);
void
check_GRC_right(double *);


int
main()
{
  // this stores the exact values of \int_0^1 x^i/sqrt(x(1-x)) dx
  static double exact_monomials[32];

  exact_monomials[0]  = 3.141592653589793;
  exact_monomials[1]  = 1.570796326794897;
  exact_monomials[2]  = 1.178097245096172;
  exact_monomials[3]  = 0.9817477042468104;
  exact_monomials[4]  = 0.8590292412159591;
  exact_monomials[5]  = 0.7731263170943632;
  exact_monomials[6]  = 0.7086991240031662;
  exact_monomials[7]  = 0.6580777580029401;
  exact_monomials[8]  = 0.6169478981277563;
  exact_monomials[9]  = 0.5826730148984365;
  exact_monomials[10] = 0.5535393641535147;
  exact_monomials[11] = 0.5283784839647186;
  exact_monomials[12] = 0.5063627137995220;
  exact_monomials[13] = 0.4868872248072327;
  exact_monomials[14] = 0.4694983953498315;
  exact_monomials[15] = 0.4538484488381705;
  exact_monomials[16] = 0.4396656848119776;
  exact_monomials[17] = 0.4267343411410371;
  exact_monomials[18] = 0.4148806094426750;
  exact_monomials[19] = 0.4039626986678677;
  exact_monomials[20] = 0.3938636312011710;
  exact_monomials[21] = 0.3844859256963813;
  exact_monomials[22] = 0.3757476092032817;
  exact_monomials[23] = 0.3675791829162538;
  exact_monomials[24] = 0.3599212832721652;
  exact_monomials[25] = 0.3527228576067219;
  exact_monomials[26] = 0.3459397257296695;
  exact_monomials[27] = 0.3395334345124534;
  exact_monomials[28] = 0.3334703374675882;
  exact_monomials[29] = 0.3277208488905608;
  exact_monomials[30] = 0.3222588347423848;
  exact_monomials[31] = 0.3170611116013786;


  initlog();
  deallog << std::setprecision(8);

  deallog << "* 1d Gauss-Chebyshev" << std::endl;
  check_quadrature<QGaussChebyshev<1>, 1>(&exact_monomials[0]);

  deallog << "* 1d Gauss-Radau-Chebyshev, left endpoint" << std::endl;
  check_quadrature<QGaussRadauChebyshev<1>, 1>(&exact_monomials[0]);

  deallog << "* 1d Gauss-Radau-Chebyshev, right endpoint" << std::endl;
  check_GRC_right(&exact_monomials[0]);

  deallog << "1d Gauss-Lobatto-Chebyshev" << std::endl;
  check_quadrature<QGaussLobattoChebyshev<1>, 2>(&exact_monomials[0]);

  return 0;
}


template <typename quadrature_type, unsigned short startn>
void
check_quadrature(double *exact_monomials)
{
  for (unsigned int n = startn; n < 18; ++n)
    {
      quadrature_type              quadrature(n);
      const std::vector<Point<1>> &points  = quadrature.get_points();
      const std::vector<double>   &weights = quadrature.get_weights();


      for (unsigned int i = 0; i < 32; ++i)
        {
          long double quadrature_int = 0;
          double      err            = 0;

          // Check the integral
          // x^i/sqrt(x(1-x))
          long double f = 1.;
          for (unsigned int x = 0; x < quadrature.size(); ++x)
            {
              f = std::pow(static_cast<long double>(points[x][0]), i * 1.0L);
              quadrature_int += f * static_cast<long double>(weights[x]);
            }
          err = std::fabs(quadrature_int - exact_monomials[i]);
          deallog << "Quadrature order " << n << ", polynomial of degree " << i
                  << ": ";

          if (err < 1.e-14)
            deallog << "exact." << std::endl;
          else
            deallog << "error " << err << std::endl;
        }
    }
}


void
check_GRC_right(double *exact_monomials)
{
  for (unsigned int n = 1; n < 18; ++n)
    {
      QGaussRadauChebyshev<1> quadrature(
        n, QGaussRadauChebyshev<1>::EndPoint::right);
      const std::vector<Point<1>> &points  = quadrature.get_points();
      const std::vector<double>   &weights = quadrature.get_weights();


      for (unsigned int i = 0; i < 32; ++i)
        {
          long double quadrature_int = 0;
          double      err            = 0;

          // Check the integral
          // x^i/sqrt(x(1-x))
          long double f = 1.;
          for (unsigned int x = 0; x < quadrature.size(); ++x)
            {
              f = std::pow(static_cast<long double>(points[x][0]), i * 1.0L);
              quadrature_int += f * static_cast<long double>(weights[x]);
            }
          err = std::fabs(quadrature_int - exact_monomials[i]);
          deallog << "Quadrature order " << n << ", polynomial of degree " << i
                  << ": ";

          if (err < 2.e-15)
            deallog << "exact." << std::endl;
          else
            deallog << "error " << err << std::endl;
        }
    }
}
