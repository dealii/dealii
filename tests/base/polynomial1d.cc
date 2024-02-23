// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// just output a lot of information about various classes implementing
// 1d-polynomials, to make sure that all changes we make to these classes
// do not change the results of these classes.

#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>

#include <iostream>

#include "../tests.h"


using namespace Polynomials;


double
scalar_product(const Polynomial<double> &p1, const Polynomial<double> &p2)
{
  unsigned int degree = (p1.degree() + p2.degree()) / 2 + 1;
  QGauss<1>    gauss(degree);

  double sum = 0.;
  for (unsigned int i = 0; i < gauss.size(); ++i)
    {
      double x  = gauss.point(i)[0];
      double P1 = p1.value(x);
      double P2 = p2.value(x);
      sum += gauss.weight(i) * P1 * P2;
    }
  return sum;
}


void
polynomial_arithmetic()
{
  std::vector<double> c1(4);
  c1[0] = 2;
  c1[1] = 3;
  c1[2] = 5;
  c1[3] = 4;
  Polynomial<double> p1(c1);

  std::vector<double> c2(3);
  c2[0] = .4;
  c2[1] = .7;
  c2[2] = -1.3;
  Polynomial<double> p2(c2);
  Monomial<double>   p3(5);
  deallog << "P1" << std::endl;
  p1.print(deallog.get_file_stream());
  deallog << "P2" << std::endl;
  p2.print(deallog.get_file_stream());
  p1 += p2;
  deallog << "P1+P2" << std::endl;
  p1.print(deallog.get_file_stream());
  p2 += p1;
  deallog << "P1+2P2" << std::endl;
  p2.print(deallog.get_file_stream());
  deallog << "P1+P2+x^5" << std::endl;
  p1 += p3;
  p1.print(deallog.get_file_stream());
  p1 *= 2.;
  deallog << "*2" << std::endl;
  p1.print(deallog.get_file_stream());
  deallog << "*P2" << std::endl;
  p2 *= p1;
  p2.print(deallog.get_file_stream());

  for (unsigned int i = 0; i < 7; ++i)
    {
      deallog << "derive" << std::endl;
      p1 = p1.derivative();
      p1.print(deallog.get_file_stream());
    }
}


int
main()
{
  initlog();

  polynomial_arithmetic();

  std::vector<Polynomial<double>> p;
  std::vector<Polynomial<double>> q;

  deallog << "Legendre" << std::endl;

  for (unsigned int i = 0; i < 12; ++i)
    {
      p.push_back(Legendre(i));
    }


  for (unsigned int i = 0; i < p.size(); ++i)
    for (unsigned int j = 0; j <= i; ++j)
      deallog << 'P' << i << " * P" << j << " = " << scalar_product(p[i], p[j])
              << std::endl;


  deallog << "LagrangeEquidistant" << std::endl;

  p.clear();
  for (unsigned int i = 0; i < 6; ++i)
    {
      p.push_back(LagrangeEquidistant(6, i));
      q.push_back(LagrangeEquidistant(6, i));
    }

  // We add 1.0001 because of bugs in
  // the ostream classes
  for (unsigned int i = 0; i < p.size(); ++i)
    for (unsigned int j = 0; j < p.size(); ++j)
      deallog << 'P' << i << "(x" << j
              << ") = " << p[i].value((double)j / p.size()) + 1.0001
              << std::endl;

  for (unsigned int i = 0; i < p.size(); ++i)
    {
      q[i].scale(.5);
      for (unsigned int j = 0; j < p.size(); ++j)
        {
          double x = (double)j / p.size();
          if (std::fabs(q[i].value(2. * x) - p[i].value(x)) > 1.e-15)
            deallog << "Polynomial " << i << ": Values q(" << 2. * x
                    << ") and p(" << x
                    << ") differ after scale: " << q[i].value(2. * x)
                    << " != " << p[i].value(x) << std::endl;
        }
      q[i].shift((double)1.);
      for (unsigned int j = 0; j < p.size(); ++j)
        {
          double x    = (double)j / p.size();
          double diff = std::fabs(q[i].value(2. * x - 1.) - p[i].value(x));
          if (diff > 1.e-13)
            deallog << "Polynomial " << i << ": Values q(" << 2. * x - 1.
                    << ") and p(" << x << ") differ by 10^"
                    << std::log(diff) / std::log(10.)
                    << " after shift: " << q[i].value(2. * x - 1.)
                    << " != " << p[i].value(x) << std::endl;
        }
    }

  deallog << "Hierarchical" << std::endl;

  p.clear();
  for (unsigned int i = 0; i < 12; ++i)
    {
      p.push_back(Hierarchical(i));
    }

  for (unsigned int i = 0; i < p.size(); ++i)
    {
      deallog << "N0(P" << i << ')' << " = " << p[i].value(0.) << std::endl;
    }

  for (unsigned int i = 0; i < p.size(); ++i)
    {
      deallog << "N1(P" << i << ')' << " = " << p[i].value(1.) << std::endl;
    }

  std::vector<double> values(p.size());

  for (unsigned int j = 2; j < p.size(); ++j)
    {
      double factor = 1.;
      for (unsigned int k = 0; k < j; ++k)
        factor /= 2. * (j - k);

      for (unsigned int i = 0; i < p.size(); ++i)
        {
          p[i].value(.5, values);
          deallog << 'N' << j << "(P" << i << ')' << " = " << values[j] * factor
                  << std::endl;
        }
    }
}
