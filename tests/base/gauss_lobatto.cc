// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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


// check points and weights for Gauss-Lobatto quadrature formula

#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"


int
main()
{
  initlog();

  for (unsigned int n = 2; n < 20; ++n)
    {
      deallog << "QGaussLobatto(" << n << ")" << std::endl;

      QGaussLobatto<1> q(n);
      for (unsigned int i = 0; i < q.size(); ++i)
        deallog << q.point(i) << ' ' << q.weight(i) << std::endl;

      // the points must be symmetrically located around 0.5
      double p = 0;
      for (unsigned int i = 0; i < q.size(); ++i)
        p += (q.point(i)[0] - 0.5);
      AssertThrow(std::fabs(p) < 1e-12, ExcInternalError());

      // the sum of weights must be one
      double w = 0;
      for (unsigned int i = 0; i < q.size(); ++i)
        w += q.weight(i);
      AssertThrow(std::fabs(w - 1) < 1e-12, ExcInternalError());
    }
}
