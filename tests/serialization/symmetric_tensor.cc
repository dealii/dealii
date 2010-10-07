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

// check serialization for SymmetricTensor<2,dim>

#include "serialization.h"
#include <base/symmetric_tensor.h>


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
  std::ofstream logfile("symmetric_tensor/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
