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

// Tests the computation of the second derivatives of a function using
// nested reverse-forward mode AD. The Sacado::Rad::ADvar class is used to
// compute the first derivatives while the Sacado::Fad::DFad class, which
// utilizes dynamic memory allocation for the number of derivative components,
// is used for the second derivative calculations.  The difference between this
// test and sacado/basic_01b.cc is that we test the vector-mode capability of
// Sacado::Rad::ADvar.
//
// A related examples that are shipped with Trilinos can be found at
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/trad_dfad_example.cpp
// https://github.com/trilinos/Trilinos/blob/master/packages/sacado/example/tradvec_example.cpp


#include <Sacado.hpp>
#include <Sacado_trad.hpp>

#include "../tests.h"

// The function to differentiate
template <typename NumberType, typename NumberType2>
NumberType
f(const NumberType &x, const NumberType &y, const NumberType2 &z)
{
  return z * (x * x * x + z * y * y + 0.5 * x * y * y);
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
  return x * x * y * y * z;
}

// The analytic derivative of the functions with respect to x and y
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
  dh_dx = 2 * x * y * y * z;
  dh_dy = 2 * x * x * y * z;
}

// The analytic second derivatives of the functions with respect to x and y
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
void
d2g(const double &x,
    const double &y,
    const double &z,
    double       &d2g_dx_dx,
    double       &d2g_dy_dy,
    double       &d2g_dy_dx)
{
  d2g_dx_dx = -z * z * std::sin(x * z) * std::cos(y / z);
  d2g_dy_dx = -std::cos(x * z) * std::sin(y / z);
  d2g_dy_dy = -(1.0 / (z * z)) * std::sin(x * z) * std::cos(y / z);
}
void
d2h(const double &x,
    const double &y,
    const double &z,
    double       &d2h_dx_dx,
    double       &d2h_dy_dy,
    double       &d2h_dy_dx)
{
  d2h_dx_dx = 2 * y * y * z;
  d2h_dy_dx = 4 * x * y * z;
  d2h_dy_dy = 2 * x * x * z;
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
  Sacado::Rad::ADvar<Sacado::Fad::DFad<double>> x_ad(
    Sacado::Fad::DFad<double>(num_deriv, 0, x));
  Sacado::Rad::ADvar<Sacado::Fad::DFad<double>> y_ad(
    Sacado::Fad::DFad<double>(num_deriv, 1, y));
  // FAD objects: Passive variables
  const Sacado::Rad::ADvar<Sacado::Fad::DFad<double>> z_ad(z);

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

  // Compute second derivative analytically
  double d2f_dx_dx = 0.0, d2f_dy_dy = 0.0, d2f_dy_dx = 0.0;
  d2f(x, y, z, d2f_dx_dx, d2f_dy_dy, d2f_dy_dx);
  double d2g_dx_dx = 0.0, d2g_dy_dy = 0.0, d2g_dy_dx = 0.0;
  d2g(x, y, z, d2g_dx_dx, d2g_dy_dy, d2g_dy_dx);
  double d2h_dx_dx = 0.0, d2h_dy_dy = 0.0, d2h_dy_dx = 0.0;
  d2h(x, y, z, d2h_dx_dx, d2h_dy_dy, d2h_dy_dx);

  // Compute function values
  // We specifically choose to do all of these computations
  // before computing gradients, because this mixes the operations
  // performed with each independent variables to produce each
  // dependent variable
  Sacado::Rad::ADvar<Sacado::Fad::DFad<double>> f_rfad =
    ::f(x_ad, y_ad, z_ad); // Cannot be const
  Sacado::Rad::ADvar<Sacado::Fad::DFad<double>> h_rfad =
    ::h(x_ad, y_ad, z_ad); // Cannot be const <----- Before g_rad
  Sacado::Rad::ADvar<Sacado::Fad::DFad<double>> g_rfad =
    ::g(x_ad, y_ad, z_ad); // Cannot be const
  deallog << "f_rfad: " << f_rfad.val() << std::endl;
  deallog << "g_rfad: " << g_rfad.val() << std::endl;
  deallog << "h_rfad: " << h_rfad.val() << std::endl;

  // Partial derivative accumulation terms
  Sacado::Fad::DFad<double> d_dx_rad_acc = 0.0;
  Sacado::Fad::DFad<double> d_dy_rad_acc = 0.0;

  // Configure the AD number to perform gradient computations
  // related to the dependent function "f"
  Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>::Outvar_Gradcomp(f_rfad);
  // Extract value and derivatives
  const double                    f_ad         = f_rfad.val().val(); // f
  const Sacado::Fad::DFad<double> df_dx_fad    = x_ad.adj();         // df/dx
  const Sacado::Fad::DFad<double> df_dy_fad    = y_ad.adj();         // df/dy
  const double                    df_dx_ad     = df_dx_fad.val();    // df/dx
  const double                    df_dy_ad     = df_dy_fad.val();    // df/dy
  const double                    d2f_dx_dx_ad = x_ad.adj().dx(0); // d^2f/dx^2
  const double                    d2f_dy_dx_ad = x_ad.adj().dx(1); // d^2f/dy_dx
  const double                    d2f_dx_dy_ad = y_ad.adj().dx(0); // d^2f/dx_dy
  const double                    d2f_dy_dy_ad = y_ad.adj().dx(1); // d^2f/dy^2

  std::cout << "df_dx: " << df_dx << "  df_dx_ad: " << df_dx_ad << std::endl;
  std::cout << "df_dy: " << df_dy << "  df_dy_ad: " << df_dy_ad << std::endl;

  // Configure the AD number to perform gradient computations
  // related to the dependent function "g"
  d_dx_rad_acc += df_dx_fad;
  d_dy_rad_acc += df_dy_fad;
  Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>::Outvar_Gradcomp(g_rfad);
  // Extract value and derivatives
  const double                    g_ad = g_rfad.val().val(); // g
  const Sacado::Fad::DFad<double> dg_dx_fad =
    x_ad.adj() -
    d_dx_rad_acc; // dg/dx ; Note: Accumulation of partial derivatives
  const Sacado::Fad::DFad<double> dg_dy_fad =
    y_ad.adj() -
    d_dy_rad_acc; // dg/dy ; Note: Accumulation of partial derivatives
  const double dg_dx_ad     = dg_dx_fad.val(); // dg/dx
  const double dg_dy_ad     = dg_dy_fad.val(); // dg/dy
  const double d2g_dx_dx_ad = dg_dx_fad.dx(0); // d^2g/dx^2
  const double d2g_dy_dx_ad = dg_dx_fad.dx(1); // d^2g/dy_dx
  const double d2g_dx_dy_ad = dg_dy_fad.dx(0); // d^2g/dx_dy
  const double d2g_dy_dy_ad = dg_dy_fad.dx(1); // d^2g/dy^2

  std::cout << "dg_dx: " << dg_dx << "  dg_dx_ad: " << dg_dx_ad << std::endl;
  std::cout << "dg_dy: " << dg_dy << "  dg_dy_ad: " << dg_dy_ad << std::endl;

  // Configure the AD number to perform gradient computations
  // related to the dependent function "h"
  d_dx_rad_acc += dg_dx_fad;
  d_dy_rad_acc += dg_dy_fad;
  Sacado::Rad::ADvar<Sacado::Fad::DFad<double>>::Outvar_Gradcomp(h_rfad);
  // Extract value and derivatives
  const double                    h_ad = h_rfad.val().val(); // h
  const Sacado::Fad::DFad<double> dh_dx_fad =
    x_ad.adj() -
    d_dx_rad_acc; // dh/dx ; Note: Accumulation of partial derivatives
  const Sacado::Fad::DFad<double> dh_dy_fad =
    y_ad.adj() -
    d_dy_rad_acc; // dh/dy ; Note: Accumulation of partial derivatives
  const double dh_dx_ad     = dh_dx_fad.val(); // dg/dx
  const double dh_dy_ad     = dh_dy_fad.val(); // dg/dy
  const double d2h_dx_dx_ad = dh_dx_fad.dx(0); // d^2h/dx^2
  const double d2h_dy_dx_ad = dh_dx_fad.dx(1); // d^2h/dy_dx
  const double d2h_dx_dy_ad = dh_dy_fad.dx(0); // d^2h/dx_dy
  const double d2h_dy_dy_ad = dh_dy_fad.dx(1); // d^2h/dy^2
  // Observation: The accumulation of the adjoints appears to be related to
  // the order in which ::Outvar_Gradcomp is called (i.e. which dependent
  // variables the adjoints are computed for), rather than the order in
  // which the functions themselves are evaluated.

  std::cout << "dh_dx: " << dh_dx << "  dh_dx_ad: " << dh_dx_ad << std::endl;
  std::cout << "dh_dy: " << dh_dy << "  dh_dy_ad: " << dh_dy_ad << std::endl;

  const double tol = 1.0e-12;
  Assert(std::fabs(f - f_ad) < tol,
         ExcMessage("Computation incorrect: Value of f"));
  Assert(std::fabs(df_dx - df_dx_ad) < tol && std::fabs(df_dy - df_dy_ad) < tol,
         ExcMessage("Computation incorrect: First derivative of f"));
  Assert(std::fabs(d2f_dx_dx - d2f_dx_dx_ad) < tol &&
           std::fabs(d2f_dy_dy - d2f_dy_dy_ad) < tol &&
           std::fabs(d2f_dy_dx - d2f_dy_dx_ad) < tol,
         ExcMessage("Computation incorrect: Second derivative of f"));
  Assert(std::fabs(g - g_ad) < tol,
         ExcMessage("Computation incorrect: Value of g"));
  Assert(std::fabs(dg_dx - dg_dx_ad) < tol && std::fabs(dg_dy - dg_dy_ad) < tol,
         ExcMessage("Computation incorrect: First derivative of g"));
  Assert(std::fabs(d2g_dx_dx - d2g_dx_dx_ad) < tol &&
           std::fabs(d2g_dy_dy - d2g_dy_dy_ad) < tol &&
           std::fabs(d2g_dy_dx - d2g_dy_dx_ad) < tol,
         ExcMessage("Computation incorrect: Second derivative of g"));
  Assert(std::fabs(h - h_ad) < tol,
         ExcMessage("Computation incorrect: Value of h"));
  Assert(std::fabs(dh_dx - dh_dx_ad) < tol && std::fabs(dh_dy - dh_dy_ad) < tol,
         ExcMessage("Computation incorrect: First derivative of h"));
  Assert(std::fabs(d2h_dx_dx - d2h_dx_dx_ad) < tol &&
           std::fabs(d2h_dy_dy - d2h_dy_dy_ad) < tol &&
           std::fabs(d2h_dy_dx - d2h_dy_dx_ad) < tol,
         ExcMessage("Computation incorrect: Second derivative of h"));

  deallog << "OK" << std::endl;
}
