// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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


// test parallel::transform

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/parallel.h>
#include <deal.II/lac/vector.h>
#include <boost/lambda/lambda.hpp>



int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const unsigned int N=10000;
  Vector<double> x(N), y(N), z(N);

  for (unsigned int i=0; i<N; ++i)
    {
      x(i) = 2.*i;
      y(i) = -1.*i;
    }

  // set z=x+2y, which happens to be zero
  parallel::transform (x.begin(), x.end(),
                       y.begin(),
                       z.begin(),
                       (boost::lambda::_1 + 2*boost::lambda::_2),
                       10);

  Assert (z.l2_norm() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}
