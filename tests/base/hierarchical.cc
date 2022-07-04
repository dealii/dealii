// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


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
