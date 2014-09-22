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


// check Tensor<0,dim>

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

template <typename U, typename V>
void compare (const U &u, const V &v)
{
  Assert (static_cast<double>(u)==static_cast<double>(v), ExcInternalError());
}



int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  typedef Tensor<0,1> T;
  T t1(13.), t2(42);

  compare (T(), 0.);
  compare (T(13.), 13.);
  compare (T(t1), 13.);
  compare (static_cast<double>(t1), 13.);
  compare (static_cast<double &>(t1), 13.);
  compare ((T() = t1), 13.);
  compare ((T() = 13.), 13.);
  compare ((t1==t1), true);
  compare ((t1==t2), false);
  compare ((t1!=t2), true);
  compare ((t1!=t1), false);
  compare ((T() += t1), t1);
  compare ((T() -= t1), -t1);
  compare (T(13.)*=3., 39.);
  compare (T(39)/=3., 13.);
  compare ((t1*t2), 13*42);
  compare ((t1+t2), 13+42);
  compare ((t1-t2), 13-42);
  compare (-t1, -13.);
  compare (T(-12).norm(), 12.);
  compare (T(-12).norm_square(), 12*12.);

  t1.clear();
  compare (t1, 0.);

  deallog << "OK" << std::endl;
}
