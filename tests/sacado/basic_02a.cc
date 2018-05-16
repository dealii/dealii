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

// Tests the computation of the second derivatives of a function using
// nested reverse-forward mode AD. The Sacado::Rad::ADvar class is used to
// compute the first derivatives while the Sacado::Fad::DFad class, which utilizes
// dynamic memory allocation for the number of derivative components, is
// used for the second derivative calculations.
//
// A related example that is shipped with Trilinos can be found at
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/trad_dfad_example.cpp


#include "../tests.h"

#include <Sacado.hpp>
#include <Sacado_trad.hpp>

// The function to differentiate
template <typename NumberType, typename NumberType2>
NumberType
f(const NumberType &x, const NumberType &y, const NumberType2 &z)
{
  return z*(x*x*x + z*y*y + 0.5*x*y*y);
}

// The analytic derivative of f(x,y,z) with respect to x and y
void
df(const double &x, const double &y, const double &z,
   double &df_dx, double &df_dy)
{
  df_dx = z*(3.0*x*x + 0.5*y*y);
  df_dy = z*(2.0*z*y + x*y);
}

// The analytic second derivatives of f(x,y,z) with respect to x and y
void
d2f(const double &x, const double &y, const double &z,
    double &d2f_dx_dx, double &d2f_dy_dy,
    double &d2f_dy_dx)
{
  d2f_dx_dx = z*(6.0*x);
  d2f_dy_dx = z*y;
  d2f_dy_dy = z*(2.0*z + x);
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
  Sacado::Rad::ADvar< Sacado::Fad::DFad<double> > x_ad(Sacado::Fad::DFad<double>(num_deriv, 0, x));
  Sacado::Rad::ADvar< Sacado::Fad::DFad<double> > y_ad(Sacado::Fad::DFad<double>(num_deriv, 1, y));
  // FAD objects: Passive variables
  const Sacado::Rad::ADvar< Sacado::Fad::DFad<double> > z_ad(z);

  deallog << "x_ad: " << x_ad.val() << std::endl;
  deallog << "y_ad: " << y_ad.val() << std::endl;
  deallog << "z_ad: " << z_ad.val() << std::endl;

  // Compute function
  const double f = ::f(x, y, z);

  // Compute derivative analytically
  double df_dx = 0.0, df_dy = 0.0;
  df(x, y, z, df_dx, df_dy);

  // Compute second derivative analytically
  double d2f_dx_dx = 0.0, d2f_dy_dy = 0.0, d2f_dy_dx = 0.0;
  d2f(x, y, z, d2f_dx_dx, d2f_dy_dy, d2f_dy_dx);

  // Compute function and derivative with AD
  const Sacado::Rad::ADvar< Sacado::Fad::DFad<double> > f_rfad = ::f(x_ad, y_ad, z_ad);
  Sacado::Rad::ADvar< Sacado::Fad::DFad<double> >::Gradcomp();

  deallog << "f_rad: " << f_rfad.val() << std::endl;

  // Extract value and derivatives
  const double f_ad = f_rfad.val().val();        // f
  const double df_dx_ad = x_ad.adj().val();      // df/dx
  const double df_dy_ad = y_ad.adj().val();      // df/dy
  const double d2f_dx_dx_ad = x_ad.adj().dx(0);  // d^2f/dx^2
  const double d2f_dy_dx_ad = x_ad.adj().dx(1);  // d^2f/dy_dx
  const double d2f_dx_dy_ad = y_ad.adj().dx(0);  // d^2f/dx_dy
  const double d2f_dy_dy_ad = y_ad.adj().dx(1);  // d^2f/dy^2

  const double tol = 1.0e-14;
  Assert(std::fabs(f - f_ad) < tol,
         ExcMessage("Computation incorrect: Value"));
  Assert(std::fabs(df_dx - df_dx_ad) < tol &&
         std::fabs(df_dy - df_dy_ad) < tol,
         ExcMessage("Computation incorrect: First derivative"));
  Assert(std::fabs(d2f_dx_dx - d2f_dx_dx_ad) < tol &&
         std::fabs(d2f_dy_dy - d2f_dy_dy_ad) < tol &&
         std::fabs(d2f_dy_dx - d2f_dy_dx_ad) < tol,
         ExcMessage("Computation incorrect: Second derivative"));

  deallog << "OK" << std::endl;
}
