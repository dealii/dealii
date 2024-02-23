// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Output function values and derivatives for vector valued
// polynomials returning Tensor<1,dim> for their values.

// classes: PolynomialsBDM, PolynomialsRaviartThomas

#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/polynomials_raviart_thomas.h>

#include <vector>

#include "../tests.h"

template <int dim, class PolynomialType>
void
check_point(const Point<dim> &x, const PolynomialType &p)
{
  const unsigned int          n = p.n();
  std::vector<Tensor<1, dim>> values(n);
  std::vector<Tensor<2, dim>> gradients(n);
  std::vector<Tensor<3, dim>> seconds(0);
  std::vector<Tensor<4, dim>> thirds(0);
  std::vector<Tensor<5, dim>> fourths(0);

  p.evaluate(x, values, gradients, seconds, thirds, fourths);

  deallog << "Point " << x << std::endl;
  for (unsigned int i = 0; i < n; ++i)
    {
      deallog << "p[" << i << "] value ";
      for (unsigned int d = 0; d < dim; ++d)
        deallog << (int)(values[i][d] + .5) << ' ';
      deallog << " first ";
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          deallog << (int)(gradients[i][d1][d2] + .5) << ' ';
      deallog << std::endl;
    }
}


template <int dim>
void
check_bdm()
{
  Point<dim> x;

  PolynomialsBDM<dim> p1(1);
  PolynomialsBDM<dim> p2(2);
  PolynomialsBDM<dim> p3(3);
  PolynomialsBDM<dim> p4(4);

  x[0] = 2.;
  x[1] = 3.;
  if (dim > 2)
    x[2] = 4;

  check_point(x, p1);
  check_point(x, p2);
  check_point(x, p3);
  check_point(x, p4);
}

template <int dim>
void
check_rt()
{
  Point<dim> x;

  PolynomialsRaviartThomas<dim> p0(0);
  PolynomialsRaviartThomas<dim> p1(1);
  PolynomialsRaviartThomas<dim> p2(2);
  PolynomialsRaviartThomas<dim> p3(3);

  x[0] = 2.;
  x[1] = 3.;
  if (dim > 2)
    x[2] = 4;

  check_point(x, p0);
  check_point(x, p1);
  check_point(x, p2);
  check_point(x, p3);
}

int
main()
{
  initlog();
  deallog << std::setprecision(0);

  deallog.push("BDM");
  check_bdm<2>();
  check_bdm<3>();
  deallog.pop();
  deallog.push("RT");
  check_rt<2>();
  check_rt<3>();
  deallog.pop();
}
