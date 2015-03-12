// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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


// Same as tensor.cc, but uses tensors based on floats instead of doubles

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  float a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};
  float b[3]    = {25,31,37};

  const unsigned int dim=3;
  
  // rank 2
  {
    Tensor<2,dim,float> t(a);
    Tensor<2,dim,double> dt(t), dt2;
    dt2 = t;
    Assert (dt2 == dt, ExcInternalError());
    Assert (dt == dt2, ExcInternalError());
  
    Tensor<2,dim,float> ft(dt), ft2;
    ft2 = dt;
    Assert (ft2 == ft, ExcInternalError());
    Assert (ft == ft2, ExcInternalError());
  
    Tensor<2,dim,std::complex<double> > ct(dt), ct2;
    ct2 = dt;
    Assert (ct2 == ct, ExcInternalError());
    Assert (ct == ct2, ExcInternalError());
  }
  
  // rank 1
  {
    Tensor<1,dim,float> t(b);
    Tensor<1,dim,double> dt(t), dt2;
    dt2 = t;
    Assert (dt2 == dt, ExcInternalError());
    Assert (dt == dt2, ExcInternalError());
  
    Tensor<1,dim,float> ft(dt),ft2;
    ft2 = dt;
    Assert (ft2 == ft, ExcInternalError());
    Assert (ft == ft2, ExcInternalError());
  
    Tensor<1,dim,std::complex<double> > ct(dt), ct2;
    ct2 = dt;
    Assert (ct2 == ct, ExcInternalError());
    Assert (ct == ct2, ExcInternalError());
  }
  
  deallog << "OK." << std::endl;
}
