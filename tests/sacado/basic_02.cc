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

// Tests the computation of the second derivatives of a function using
// nested forward mode AD. The Sacado::Fad::DFad class, which utilizes
// dynamic memory allocation for the number of derivative components, is
// used for both tiers of the derivative calculations.
//
// A related example that is shipped with Trilinos can be found at
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/dfad_dfad_example.cpp


#include <Sacado.hpp>

#include "../tests.h"

// The function to differentiate
template <typename NumberType, typename NumberType2>
NumberType
f(const NumberType &x, const NumberType &y, const NumberType2 &z)
{
  return z * (x * x * x + z * y * y + 0.5 * x * y * y);
}

// The analytic derivative of f(x,y,z) with respect to x and y
void
df(const double &x,
   const double &y,
   const double &z,
   double       &df_dx,
   double       &df_dy)
{
  df_dx = z * (3.0 * x * x + 0.5 * y * y);
  df_dy = z * (2.0 * z * y + x * y);
}

// The analytic second derivatives of f(x,y,z) with respect to x and y
void
d2f(const double &x,
    const double &y,
    const double &z,
    double       &d2f_dx_dx,
    double       &d2f_dy_dy,
    double       &d2f_dy_dx)
{
  d2f_dx_dx = z * (6.0 * x);
  d2f_dy_dx = z * y;
  d2f_dy_dy = z * (2.0 * z + x);
}

int
main()
{
  initlog();

  // Values of function arguments
  const double x = -3.0;
  const double y = 2.0;
  const double z = 7.0;

  // Number of independent variables
  const int num_deriv = 2;

  // FAD objects: Independent variables
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> x_ad(num_deriv, 0, x);
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> y_ad(num_deriv, 1, y);
  // FAD objects: Passive variables
  const Sacado::Fad::DFad<double> z_ad(z);

  // Initialize the internal data of variables from which
  // second derivatives will be computed
  x_ad.val() = Sacado::Fad::DFad<double>(num_deriv, 0, x);
  y_ad.val() = Sacado::Fad::DFad<double>(num_deriv, 1, y);

  deallog << "x_ad: " << x_ad << std::endl;
  deallog << "y_ad: " << y_ad << std::endl;
  deallog << "z_ad: " << z_ad << std::endl;

  // Compute function
  const double f = ::f(x, y, z);

  // Compute derivative analytically
  double df_dx = 0.0, df_dy = 0.0;
  df(x, y, z, df_dx, df_dy);

  // Compute second derivative analytically
  double d2f_dx_dx = 0.0, d2f_dy_dy = 0.0, d2f_dy_dx = 0.0;
  d2f(x, y, z, d2f_dx_dx, d2f_dy_dy, d2f_dy_dx);

  // Compute function and derivative with AD
  const Sacado::Fad::DFad<Sacado::Fad::DFad<double>> f_fad =
    ::f(x_ad, y_ad, z_ad);

  deallog << "f_fad: " << f_fad << std::endl;

  // Extract value and derivatives
  const double f_ad         = f_fad.val().val(); // f
  const double df_dx_ad     = f_fad.dx(0).val(); // df/dx
  const double df_dy_ad     = f_fad.dx(1).val(); // df/dy
  const double d2f_dx_dx_ad = f_fad.dx(0).dx(0); // d^2f/dx^2
  const double d2f_dy_dx_ad = f_fad.dx(0).dx(1); // d^2f/dy_dx
  const double d2f_dx_dy_ad = f_fad.dx(1).dx(0); // d^2f/dx_dy
  const double d2f_dy_dy_ad = f_fad.dx(1).dx(1); // d^2f/dy^2

  const double tol = 1.0e-14;
  Assert(std::fabs(f - f_ad) < tol, ExcMessage("Computation incorrect: Value"));
  Assert(std::fabs(df_dx - df_dx_ad) < tol && std::fabs(df_dy - df_dy_ad) < tol,
         ExcMessage("Computation incorrect: First derivative"));
  Assert(std::fabs(d2f_dx_dx - d2f_dx_dx_ad) < tol &&
           std::fabs(d2f_dy_dy - d2f_dy_dy_ad) < tol &&
           std::fabs(d2f_dy_dx - d2f_dy_dx_ad) < tol,
         ExcMessage("Computation incorrect: Second derivative"));

  deallog << "OK" << std::endl;
}
