// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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


// Tests operations similar to polynomial_lagrange_order when the polynomial
// is modified

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/logstream.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>


using namespace Polynomials;


void check_scale (const std::vector<Polynomial<double> > &p)
{
  deallog << "Scale operation";
  for (unsigned int i=0; i<p.size(); ++i)
    {
      Polynomial<double> q = p[i];
      double x = (double)Testing::rand()/RAND_MAX;
      double factor = 5.*(double)Testing::rand()/RAND_MAX;
      q.scale (factor);
      double value1 = p[i].value (factor * x);
      double value2 = q   .value (x);
      if (std::fabs(value1-value2) > std::max(1e-13,1e-13*std::fabs(value1)))
        deallog << "Error scale at x=" << x << ": p(t)=" << value1
                << ", q(x)=" << value2 << std::endl;
      deallog << ".";
    }
  deallog << std::endl;
}



void check_shift (const std::vector<Polynomial<double> > &p)
{
  // shift does not work for too high orders
  if (p.size() > 30)
    return;

  deallog << "Shift operation";
  for (unsigned int i=0; i<p.size(); ++i)
    {
      Polynomial<double> q = p[i];
      double x = (double)Testing::rand()/RAND_MAX;
      double a = 10.*(-1.+2.*(double)Testing::rand()/RAND_MAX);
      q.shift (a);
      double value1 = p[i].value (x+a);
      double value2 = q   .value (x);
      if (std::fabs(value1-value2) > std::max(1e-13,1e-13*std::fabs(value1)))
        deallog << "Error shift at x=" << x << ": p(t)=" << value1
                << ", q(x)=" << value2 << std::endl;
      deallog << ".";
    }
  deallog << std::endl;
}



void check_mult_scalar (const std::vector<Polynomial<double> > &p)
{
  deallog << "Multiply by scalar";
  for (unsigned int i=0; i<p.size(); ++i)
    {
      Polynomial<double> q = p[i];
      double x = (double)Testing::rand()/RAND_MAX;
      double a = (double)Testing::rand()/RAND_MAX;
      q *= a;
      double value1 = p[i].value (x) * a;
      double value2 = q   .value (x);
      if (std::fabs(value1-value2) > std::max(1e-13,1e-13*std::fabs(value1)))
        deallog << "Error multiply at x=" << x << ": a*p(x)=" << value1
                << ", q(x)=" << value2 << std::endl;
      deallog << ".";
    }
  deallog << std::endl;
}



void check_mult (const std::vector<Polynomial<double> > &p)
{
  deallog << "Multiply by polynomial";
  for (unsigned int i=0; i<p.size(); ++i)
    {
      for (unsigned int j=0; j<p.size(); ++j)
        {
          Polynomial<double> q = p[i];
          q *= p[j];
          double x = (double)Testing::rand()/RAND_MAX;
          double value1 = p[i].value (x) * p[j].value(x);
          double value2 = q   .value (x);
          if (std::fabs(value1-value2) > std::max(1e-13,1e-13*std::fabs(value1)))
            deallog << "Error multiply at x=" << x << ": p_1(x)*p_2(x)=" << value1
                    << ", q(x)=" << value2 << std::endl;
        }
      deallog << ".";
    }
  deallog << std::endl;
}



void check_expand (const std::vector<Polynomial<double> > &p)
{
  if (p.size() > 10)
    return;
  // this checks whether the Lagrange product
  // form and the usual form with factors for
  // different powers does the same
  // thing. Realize this by adding the
  // polynomial '0' to the current
  // polynomial. That destroys the product form
  // and triggers the usual form. do not do this
  // for higher order because then the standard
  // form is unstable
  deallog << "Expansion operation";
  Monomial<double> zero(0, 0.);
  for (unsigned int i=0; i<p.size(); ++i)
    {
      Polynomial<double> q = p[i];
      double x = (double)Testing::rand()/RAND_MAX;
      q += zero;
      double value1 = p[i].value (x);
      double value2 = q   .value (x);
      if (std::fabs(value1-value2) > std::max(1e-10,1e-10*std::fabs(value1)))
        deallog << "Error expansion at x=" << x << ": p(x)=" << value1
                << ", q(x)=" << value2 << std::endl;
      deallog << ".";
    }
  deallog << std::endl;
}



void check_mult_expand (const std::vector<Polynomial<double> > &p)
{
  if (p.size() > 6)
    return;
  deallog << "Multiply by polynomial expanding";
  Monomial<double> zero(0, 0.);
  for (unsigned int i=0; i<p.size(); ++i)
    {
      for (unsigned int j=0; j<p.size(); ++j)
        {
          Polynomial<double> q = p[i];
          q += zero;
          q *= p[j];
          double x = (double)Testing::rand()/RAND_MAX;
          double value1 = p[i].value (x) * p[j].value(x);
          double value2 = q   .value (x);
          if (std::fabs(value1-value2) > std::max(1e-10,1e-10*std::fabs(value1)))
            deallog << "Error multiply at x=" << x
                    << ": p_"<<i<< "(x)*p_"<<j<<"(x)=" << value1
                    << ", q(x)=" << value2 << std::endl;
        }
      deallog << ".";
    }
  deallog << std::endl;
}



void
check_lge (unsigned int n)
{
  deallog << "Points: " << n+1 << std::endl;
  std::vector<Polynomial<double> > p = LagrangeEquidistant::generate_complete_basis(n);
  check_scale (p);
  check_shift (p);
  check_mult_scalar(p);
  check_mult  (p);
  check_expand(p);
  check_mult_expand(p);
  deallog << std::endl;
}



void
check_poly (const Quadrature<1> &q)
{
  deallog << "Points: " << q.size() << std::endl;
  std::vector<Polynomial<double> > p = generate_complete_Lagrange_basis(q.get_points());
  check_scale (p);
  check_shift (p);
  check_mult_scalar(p);
  check_mult  (p);
  check_expand(p);
  check_mult_expand(p);
  deallog << std::endl;
}



int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("LagrangeEquidistant");
  for (unsigned i=1; i<10; i+=2)
    check_lge (i);
  deallog.pop();
  deallog << std::endl;

  // Lagrange elements on GL points have good
  // conditioning, so test to some very high
  // orders
  deallog.push("GaussLobatto");
  for (unsigned i=1; i<40; i+=3)
    check_poly (QGaussLobatto<1>(i+1));
  deallog.pop();
}
