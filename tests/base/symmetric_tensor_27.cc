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


// test scalar_product between tensors and symmetric tensors

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  // initialize symmetric and non-symmetric tensors. in the former case, we
  // only need to initialize one half
  SymmetricTensor<2,dim> s;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      s[i][j] = (1.+(i+1)*(j*2));

  Tensor<2,dim> t;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      t[i][j] = (1.+(i+1)*(j*2));

  deallog << scalar_product(s,s) << std::endl;
  Assert (scalar_product(s,s) == s*s, ExcInternalError());
  
  deallog << scalar_product(s,t) << std::endl;
  Assert (scalar_product(s,t) == s*symmetrize(t), ExcInternalError());
  Assert (scalar_product(s,t) == symmetrize(t)*s, ExcInternalError());
  
  deallog << scalar_product(t,s) << std::endl;
  Assert (scalar_product(t,s) == s*symmetrize(t), ExcInternalError());
  Assert (scalar_product(t,s) == symmetrize(t)*s, ExcInternalError());
}




int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
