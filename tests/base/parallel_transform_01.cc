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
  Vector<double> x(N), y(N);

  for (unsigned int i=0; i<N; ++i)
    x(i) = i;

  // set y=2*x
  parallel::transform (x.begin(), x.end(), y.begin(),
                       (2*boost::lambda::_1),
                       10);

  // compute y=0 from the previous result
  y -= x;
  y -= x;

  Assert (y.l2_norm() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}
