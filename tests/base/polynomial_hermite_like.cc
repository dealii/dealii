// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/polynomial.h>

#include "../tests.h"

using namespace Polynomials;

void
plot(const std::vector<Polynomial<double>> &polynomials)
{
  LogStream::Prefix("plot");
  const unsigned int n = 8;
  for (unsigned int i = 0; i <= n; ++i)
    {
      const double x = 1. * i / n;
      deallog << x;
      for (unsigned int p = 0; p < polynomials.size(); ++p)
        deallog << '\t' << polynomials[p].value(x);
      deallog << std::endl;
    }
}

void
interpolation_conditions(const std::vector<Polynomial<double>> &polynomials)
{
  std::vector<double> values(2);
  for (unsigned int i = 0; i < polynomials.size(); ++i)
    {
      polynomials[i].value(0., values);
      deallog << i << "\t 0:\t" << values[0] << '\t' << values[1] << std::endl;
      polynomials[i].value(1., values);
      deallog << i << "\t 1:\t" << values[0] << '\t' << values[1] << std::endl;
    }
}


int
main()
{
  initlog();

  deallog << "degree 1" << std::endl;
  plot(HermiteLikeInterpolation::generate_complete_basis(1));
  deallog << std::endl << "degree 2" << std::endl;
  plot(HermiteLikeInterpolation::generate_complete_basis(2));
  deallog << std::endl << "degree 3" << std::endl;
  interpolation_conditions(
    HermiteLikeInterpolation::generate_complete_basis(3));
  plot(HermiteLikeInterpolation::generate_complete_basis(3));
  deallog << std::endl << "degree 4" << std::endl;
  interpolation_conditions(
    HermiteLikeInterpolation::generate_complete_basis(4));
  plot(HermiteLikeInterpolation::generate_complete_basis(4));
  deallog << std::endl << "degree 6" << std::endl;
  interpolation_conditions(
    HermiteLikeInterpolation::generate_complete_basis(6));
  plot(HermiteLikeInterpolation::generate_complete_basis(6));
  deallog << std::endl << "degree 9" << std::endl;
  interpolation_conditions(
    HermiteLikeInterpolation::generate_complete_basis(9));
  plot(HermiteLikeInterpolation::generate_complete_basis(9));
}
