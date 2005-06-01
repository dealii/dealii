//----------------------------  symmetric_tensor_16.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  symmetric_tensor_16.cc  ---------------------------

// compute the deviator tensor as stated in the documentation of outer_product

#include "../tests.h"
#include <base/symmetric_tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  const SymmetricTensor<4,dim>
    T = (identity_tensor<dim>()
         - 1./dim * outer_product(unit_symmetric_tensor<dim>(),
                                  unit_symmetric_tensor<dim>()));

  Assert ((T-deviator_tensor<dim>()).norm()
          <= 1e-15*T.norm(), ExcInternalError());
}

  


int main ()
{
  std::ofstream logfile("symmetric_tensor_16.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
