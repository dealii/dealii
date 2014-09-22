// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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


// check points and weights for Gauss-Lobatto quadrature formula

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <cmath>


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  for (unsigned int n=2; n<20; ++n)
    {
      deallog << "QGaussLobatto(" << n << ")" << std::endl;

      QGaussLobatto<1> q(n);
      for (unsigned int i=0; i<q.size(); ++i)
        deallog << q.point(i) << ' ' << q.weight(i) << std::endl;

      // the points must be symmetrically located around 0.5
      double p = 0;
      for (unsigned int i=0; i<q.size(); ++i)
        p += (q.point(i)[0] - 0.5);
      Assert (std::fabs(p) < 1e-12, ExcInternalError());

      // the sum of weights must be one
      double w = 0;
      for (unsigned int i=0; i<q.size(); ++i)
        w += q.weight(i);
      Assert (std::fabs(w-1) < 1e-12, ExcInternalError());
    }
}


