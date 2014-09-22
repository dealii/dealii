// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


// check Tensor<1,dim>::operator*(Tensor<1,dim>)

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>

template <int dim>
void test_tensor ()
{
  Tensor<1,dim> t1, t2;
  for (unsigned int i=0; i<dim; ++i)
    {
      t1[i] = i+1;
      t2[i] = 4.*i-10.;
    }
  double res = t1 * t2;
  deallog << "dim = " << dim << ": " << res << std::endl;
}

int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_tensor<1>();
  test_tensor<2>();
  test_tensor<3>();
  test_tensor<4>();
  test_tensor<7>();
  deallog << "OK" << std::endl;
}
