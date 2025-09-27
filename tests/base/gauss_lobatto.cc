// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check points and weights for Gauss-Lobatto quadrature formula

#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"


int
main()
{
  initlog();

  for (unsigned int n = 2; n < 20; ++n)
    {
      deallog << "QGaussLobatto(" << n << ')' << std::endl;

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
