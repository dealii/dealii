//----------------------------  tensor_base.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_base.cc  ---------------------------

// check serialization for Tensor<1,dim>

#include "serialization.h"
#include <base/tensor.h>


void test ()
{
  const unsigned int dim=3;

  double a1[3] = {1, 2, 3};
  Tensor<1,dim> t1(a1);

  double a2[3] = {3, 6, 9};
  Tensor<1,dim> t2(a2);

  verify (t1, t2);
}


int main ()
{
  std::ofstream logfile("tensor_base/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
