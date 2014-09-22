// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


// check serialization for SymmetricTensor<2,dim>

#include "serialization.h"
#include <deal.II/base/symmetric_tensor.h>


void test ()
{
  const unsigned int dim=3;
  const unsigned int rank=2;

  double a1[3][3] = {{1., 2., 3.},
    {2., 5., 6.},
    {3., 6., 9.}
  };
  SymmetricTensor<rank,dim> t1(a1);


  double a2[3][3] = {{10., 11., 12.},
    {11., 14., 15.},
    {12., 15., 18.}
  };
  SymmetricTensor<rank,dim> t2(a2);

  verify (t1, t2);
}


int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
