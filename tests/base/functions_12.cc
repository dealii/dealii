// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test Polynomial

#include <deal.II/base/function_lib.h>

#include "../tests.h"

template <int dim>
void
check()
{
  unsigned int n_mon = 3;

  Table<2, double> exponents(n_mon, dim);

  for (unsigned int i = 0; i < n_mon; ++i)
    for (unsigned int d = 0; d < dim; ++d)
      exponents[i][d] = i + d;

  std::vector<double> coeffs(n_mon);
  for (unsigned int i = 0; i < n_mon; ++i)
    coeffs[i] = std::pow(-1.0, static_cast<double>(i)) * (i + 1);

  Functions::Polynomial<dim> poly(exponents, coeffs);

  Point<dim> p;
  for (unsigned int d = 0; d < dim; ++d)
    p[d] = d;

  deallog << dim << "-D check" << std::endl;
  deallog << "Polynomial: ";
  for (unsigned int i = 0; i < n_mon; ++i)
    {
      deallog << coeffs[i];
      for (unsigned int d = 0; d < dim; ++d)
        deallog << " x" << d << '^' << exponents[i][d];
      if (i < n_mon - 1)
        deallog << " + ";
    }
  deallog << std::endl;
  deallog << "Point: ";
  for (unsigned int d = 0; d < dim; ++d)
    deallog << p[d] << ' ';
  deallog << std::endl;

  deallog << "Value: " << poly.value(p) << std::endl;
  deallog << " values checked" << std::endl;

  deallog << "Gradient: " << poly.gradient(p) << std::endl;
  deallog << " gradients checked" << std::endl;
  deallog << std::endl;
}

int
main()
{
  initlog();

  check<1>();
  check<2>();
  check<3>();
}
