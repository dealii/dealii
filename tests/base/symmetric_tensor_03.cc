// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// test symmetric 2x2x2x2 tensors

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  SymmetricTensor<4,2> t;
  t[0][0][0][0] = 1;
  t[1][1][1][1] = 2;
  t[0][1][0][1] = 3;

  Assert (t[0][1][0][1] == t[1][0][1][0], ExcInternalError());

  // check that if a single element is
  // accessed, its transpose element gets the
  // same value
  t[1][0][0][1] = 4;
  Assert (t[0][1][1][0] == 4, ExcInternalError());

  // make sure transposition doesn't change
  // anything
  Assert (t == transpose(t), ExcInternalError());

  // check norm of tensor
  deallog << t.norm() << std::endl;

  // make sure norm is induced by scalar
  // product
  double norm_sqr = 0;
  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<2; ++j)
      for (unsigned int k=0; k<2; ++k)
        for (unsigned int l=0; l<2; ++l)
          norm_sqr += t[i][j][k][l] * t[i][j][k][l];

  Assert (std::fabs (t.norm()*t.norm() - norm_sqr) < 1e-14,
          ExcInternalError());

  deallog << "OK" << std::endl;
}
