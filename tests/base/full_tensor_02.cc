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


// test full 3x3 tensors

#include "../tests.h"
#include <deal.II/base/tensor.h>
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

  Tensor<2,3> t;
  t[0][0] = 1;
  t[1][1] = 2;
  t[2][2] = 3;
  t[0][1] = 4;
  t[1][0] = 4;
  t[0][2] = 5;
  t[2][0] = 5;
  t[1][2] = 6;
  t[2][1] = 6;

  Assert (t[0][1] == t[1][0], ExcInternalError());

  // make sure transposition doesn't change
  // anything
  Assert (t == transpose(t), ExcInternalError());

  // check norm of tensor
  Assert (std::fabs(t.norm() - std::sqrt(1.*1+2*2+3*3+2*4*4+2*5*5+2*6*6))
          < 1e-14,
          ExcInternalError());

  deallog << "OK" << std::endl;
}
