//----------------------------  symmetric_tensor_14.cc  ---------------------------
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
//----------------------------  symmetric_tensor_14.cc  ---------------------------

// test identity_tensor

#include "../tests.h"
#include <base/symmetric_tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;
  
  SymmetricTensor<2,dim> t;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      t[i][j] = (1.+(i+1)*(j*2));

  SymmetricTensor<2,dim> x = identity_tensor<dim>() * t;
  Assert ((x-t).norm() < 1e-15*t.norm(), ExcInternalError());

  deallog << "x=" << std::endl;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      deallog << i << ' ' << j << ' ' << x[i][j] << std::endl;
}

  


int main ()
{
  std::ofstream logfile("symmetric_tensor_14.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
