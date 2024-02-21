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
// forward mode AD. The Sacado::Fad::DFad class, which utilizes dynamic
// memory allocation for the number of derivative components, is selected
// as the number type with which to perform the derivative calculations.
//
// A related example that is shipped with Trilinos can be found at
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/dfad_example.cpp


#include <Sacado.hpp>

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

  // Number of independent variables
  const int num_deriv = 2;

  // FAD objects: Independent variables
  const Sacado::Fad::DFad<double> x_ad(num_deriv, 0, x);
  const Sacado::Fad::DFad<double> y_ad(num_deriv, 1, y);
  // FAD objects: Passive variables
  const Sacado::Fad::DFad<double> z_ad(z);

  deallog << "x_ad: " << x_ad << std::endl;
  deallog << "y_ad: " << y_ad << std::endl;
  deallog << "z_ad: " << z_ad << std::endl;

  // Compute function
  const double f = ::f(x, y, z);

  // Compute derivative analytically
  double df_dx = 0.0, df_dy = 0.0;
  df(x, y, z, df_dx, df_dy);

  // Compute function and derivative with AD
  const Sacado::Fad::DFad<double> f_fad = ::f(x_ad, y_ad, z_ad);

  deallog << "f_fad: " << f_fad << std::endl;

  // Extract value and derivatives
  const double f_ad     = f_fad.val(); // f
  const double df_dx_ad = f_fad.dx(0); // df/dx
  const double df_dy_ad = f_fad.dx(1); // df/dy

  const double tol = 1.0e-14;
  Assert(std::fabs(f - f_ad) < tol, ExcMessage("Computation incorrect: Value"));
  Assert(std::fabs(df_dx - df_dx_ad) < tol && std::fabs(df_dy - df_dy_ad) < tol,
         ExcMessage("Computation incorrect: First derivative"));

  deallog << "OK" << std::endl;
}
