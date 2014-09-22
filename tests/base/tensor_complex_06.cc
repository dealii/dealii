// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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


// check Tensor<1,dim,std::complex<double>>::operator*(Tensor<1,dim,std::complex<double> >)
// These tests are simple enough to do by hand.

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>

template <int dim>
void test_tensor_01 ()
{
  Tensor<1,dim,std::complex<double> > t (false);
  for (unsigned int i=0; i<dim; ++i)
    {
      t[i] = std::complex<double> (i,1);
    }
  const Tensor<1,dim,std::complex<double> > res = 
    t * std::complex<double> (1,2);

  deallog << "dim = " << dim << ": " << res << std::endl;
}


template <int dim>
void test_tensor_02 ()
{
  Tensor<1,dim,std::complex<double> > t (false);
  for (unsigned int i=0; i<dim; ++i)
    {
      t[i] = std::complex<double> (1,i);
    }
  const Tensor<1,dim,std::complex<double> > res = 
    t * std::complex<double> (1,2);

  deallog << "dim = " << dim << ": " << res << std::endl;
}

int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << "First test" << std::endl;

  test_tensor_01<1>();
  test_tensor_01<2>();
  test_tensor_01<3>();

  deallog << "Second test" << std::endl;

  test_tensor_02<1>();
  test_tensor_02<2>();
  test_tensor_02<3>();

  deallog << "OK" << std::endl;
}
