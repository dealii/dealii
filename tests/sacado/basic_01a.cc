// ------------------------------------------------------------------------
//
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

// Tests the computation of the first derivatives of a function using
// reverse mode AD with the Sacado::Rad::ADvar class.
//
// A related example that is shipped with Trilinos can be found at
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/trad_example.cpp


#include <Sacado.hpp>
#include <Sacado_trad.hpp>

#include "../tests.h"

// The function to differentiate
template <typename NumberType>
NumberType
f(const NumberType &x, const NumberType &y, const NumberType &z)
{
  return z * (x + z * y + x * y);
}

// The analytic derivative of f(x,y,z) with respect to x and y
void
df(const double &x,
   const double &y,
   const double &z,
   double       &df_dx,
   double       &df_dy)
{
  df_dx = z * (1.0 + y);
  df_dy = z * (z + x);
}

int
main()
{
  initlog();

  // Values of function arguments
  const double x = 5.0;
  const double y = 10.0;
  const double z = 4.0;

  // RAD objects: Independent variables
  const Sacado::Rad::ADvar<double> x_ad(x);
  const Sacado::Rad::ADvar<double> y_ad(y);
  // RAD objects: Passive variables
  const Sacado::Rad::ADvar<double> z_ad(z);

  deallog << "x_ad: " << x_ad.val() << std::endl;
  deallog << "y_ad: " << y_ad.val() << std::endl;
  deallog << "z_ad: " << z_ad.val() << std::endl;

  // Compute function
  const double f = ::f(x, y, z);

  // Compute derivative analytically
  double df_dx = 0.0, df_dy = 0.0;
  df(x, y, z, df_dx, df_dy);

  // Compute function and derivative with AD
  const Sacado::Rad::ADvar<double> f_rad = ::f(x_ad, y_ad, z_ad);
  Sacado::Rad::ADvar<double>::Gradcomp();

  deallog << "f_rad: " << f_rad.val() << std::endl;

  // Extract value and derivatives
  const double f_ad     = f_rad.val(); // f
  const double df_dx_ad = x_ad.adj();  // df/dx
  const double df_dy_ad = y_ad.adj();  // df/dy

  const double tol = 1.0e-14;
  Assert(std::fabs(f - f_ad) < tol, ExcMessage("Computation incorrect: Value"));
  Assert(std::fabs(df_dx - df_dx_ad) < tol && std::fabs(df_dy - df_dy_ad) < tol,
         ExcMessage("Computation incorrect: First derivative"));

  deallog << "OK" << std::endl;
}
