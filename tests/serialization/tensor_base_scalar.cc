//----------------------------  tensor_base_scalar.cc  ---------------------------
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
//----------------------------  tensor_base_scalar.cc  ---------------------------

// check serialization for Tensor<0,dim>

#include "serialization.h"
#include "../tests.h"
#include <deal.II/base/tensor.h>


void test ()
{
  const unsigned int dim=3;

  Tensor<0,dim> t1(1.);

  Tensor<0,dim> t2(2.);

  verify (t1, t2);
}


int main ()
{
  std::ofstream logfile("tensor_base_scalar/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
