// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
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
// reverse mode AD with the Sacado::Rad::ADvar class. The difference between
// this test and sacado/basic_01a.cc is that we test the vector-mode capability
// of this number type.
//
// A related example that is shipped with Trilinos can be found at
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/tradvec_example.cpp


#include <Sacado.hpp>
#include <Sacado_trad.hpp>

#include "../tests.h"

// The functions to differentiate
template <typename NumberType>
NumberType
f(const NumberType &x, const NumberType &y, const NumberType &z)
{
  return z * (x + z * y + x * y);
}
template <typename NumberType>
NumberType
g(const NumberType &x, const NumberType &y, const NumberType &z)
{
  return std::sin(x * z) * std::cos(y / z);
}
template <typename NumberType>
NumberType
h(const NumberType &x, const NumberType &y, const NumberType &z)
{
  return x * y * z;
}

// The analytic derivative of the functions with respect to x and y
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
void
dg(const double &x,
   const double &y,
   const double &z,
   double       &dg_dx,
   double       &dg_dy)
{
  dg_dx = z * std::cos(x * z) * std::cos(y / z);
  dg_dy = -(1.0 / z) * std::sin(x * z) * std::sin(y / z);
}
void
dh(const double &x,
   const double &y,
   const double &z,
   double       &dh_dx,
   double       &dh_dy)
{
  dh_dx = y * z;
  dh_dy = x * z;
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

  // Compute functions
  const double f = ::f(x, y, z);
  const double g = ::g(x, y, z);
  const double h = ::h(x, y, z);

  // Compute derivatives analytically
  double df_dx = 0.0, df_dy = 0.0;
  df(x, y, z, df_dx, df_dy);
  double dg_dx = 0.0, dg_dy = 0.0;
  dg(x, y, z, dg_dx, dg_dy);
  double dh_dx = 0.0, dh_dy = 0.0;
  dh(x, y, z, dh_dx, dh_dy);

  // Compute function values
  // We specifically choose to do all of these computations
  // before computing gradients, because this mixes the operations
  // performed with each independent variables to produce each
  // dependent variable
  Sacado::Rad::ADvar<double> f_rad = ::f(x_ad, y_ad, z_ad); // Cannot be const
  Sacado::Rad::ADvar<double> h_rad =
    ::h(x_ad, y_ad, z_ad); // Cannot be const <----- Before g_rad
  Sacado::Rad::ADvar<double> g_rad = ::g(x_ad, y_ad, z_ad); // Cannot be const
  deallog << "f_rad: " << f_rad.val() << std::endl;
  deallog << "g_rad: " << g_rad.val() << std::endl;
  deallog << "h_rad: " << h_rad.val() << std::endl;

  // Configure the AD number to perform gradient computations
  // related to the dependent function "f"
  Sacado::Rad::ADvar<double>::Outvar_Gradcomp(f_rad);
  // Extract value and derivatives
  const double f_ad     = f_rad.val(); // f
  const double df_dx_ad = x_ad.adj();  // df/dx
  const double df_dy_ad = y_ad.adj();  // df/dy

  std::cout << "df_dx: " << df_dx << "  df_dx_ad: " << df_dx_ad << std::endl;
  std::cout << "df_dy: " << df_dy << "  df_dy_ad: " << df_dy_ad << std::endl;

  // Configure the AD number to perform gradient computations
  // related to the dependent function "g"
  Sacado::Rad::ADvar<double>::Outvar_Gradcomp(g_rad);
  // Extract value and derivatives
  const double g_ad = g_rad.val(); // g
  const double dg_dx_ad =
    (x_ad.adj() -
     df_dx_ad); // dg/dx ; Note: Accumulation of partial derivatives
  const double dg_dy_ad =
    (y_ad.adj() -
     df_dy_ad); // dg/dy ; Note: Accumulation of partial derivatives

  std::cout << "dg_dx: " << dg_dx << "  dg_dx_ad: " << dg_dx_ad << std::endl;
  std::cout << "dg_dy: " << dg_dy << "  dg_dy_ad: " << dg_dy_ad << std::endl;

  // Configure the AD number to perform gradient computations
  // related to the dependent function "h"
  Sacado::Rad::ADvar<double>::Outvar_Gradcomp(h_rad);
  // Extract value and derivatives
  const double h_ad = h_rad.val(); // h
  const double dh_dx_ad =
    (x_ad.adj() - dg_dx_ad -
     df_dx_ad); // dh/dx ; Note: Accumulation of partial derivatives
  const double dh_dy_ad =
    (y_ad.adj() - dg_dy_ad -
     df_dy_ad); // dh/dy ; Note: Accumulation of partial derivatives
  // Observation: The accumulation of the adjoints appears to be related to
  // the order in which ::Outvar_Gradcomp is called (i.e. which dependent
  // variables the adjoints are computed for), rather than the order in
  // which the functions themselves are evaluated.

  std::cout << "dh_dx: " << dh_dx << "  dh_dx_ad: " << dh_dx_ad << std::endl;
  std::cout << "dh_dy: " << dh_dy << "  dh_dy_ad: " << dh_dy_ad << std::endl;

  const double tol = 1.0e-14;
  Assert(std::fabs(f - f_ad) < tol,
         ExcMessage("Computation incorrect: Value of f"));
  Assert(std::fabs(df_dx - df_dx_ad) < tol && std::fabs(df_dy - df_dy_ad) < tol,
         ExcMessage("Computation incorrect: First derivative of f"));
  Assert(std::fabs(g - g_ad) < tol,
         ExcMessage("Computation incorrect: Value of g"));
  Assert(std::fabs(dg_dx - dg_dx_ad) < tol && std::fabs(dg_dy - dg_dy_ad) < tol,
         ExcMessage("Computation incorrect: First derivative of g"));
  Assert(std::fabs(h - h_ad) < tol,
         ExcMessage("Computation incorrect: Value of h"));
  Assert(std::fabs(dh_dx - dh_dx_ad) < tol && std::fabs(dh_dy - dh_dy_ad) < tol,
         ExcMessage("Computation incorrect: First derivative of h"));

  deallog << "OK" << std::endl;
}
