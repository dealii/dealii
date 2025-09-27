// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// generate a hierarchical basis and display the values of the shape
// functions at equidistant points. (I needed this output at one point
// in time, so why not make it a testcase -- WB)

#include <deal.II/base/polynomial.h>

#include "../tests.h"


using namespace Polynomials;


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  const std::vector<Polynomial<double>> p =
    Hierarchical::generate_complete_basis(10);

  const unsigned int div = 30;
  for (unsigned int i = 0; i <= div; ++i)
    {
      const double x = 1. * i / div;
      deallog << x << ' ';
      for (unsigned int j = 0; j < p.size(); ++j)
        deallog << p[j].value(x) << ' ';
      deallog << std::endl;
    }
}
