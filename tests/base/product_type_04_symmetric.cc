// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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


// test that the product between a symmetric tensor and an integer behaves just
// like that of the tensor and the integer-converted-to-double

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <complex>

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/symmetric_tensor.h>



int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  {
    SymmetricTensor<2,2> t;
    t[0][0] = 1.23456;
    t[0][1] = 7.87965;
    t[1][1] = 3.35792;
    AssertThrow (7*t == 7.0 * t, ExcInternalError());
    AssertThrow (t*7 == t * 7.0, ExcInternalError());
    AssertThrow (t*7 == 7*t    , ExcInternalError());
    AssertThrow ((t*7 - (t+t+t+t+t+t+t)).norm() < 1e-12, ExcInternalError());
  }

  deallog << "OK" << std::endl;
}
