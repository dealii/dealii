//----------------------------  tensor_02.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_02.cc  ---------------------------

// check Tensor::operator= (double)

#include <base/tensor.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <fstream>
#include <iostream>

int main ()
{
  std::ofstream logfile("tensor_02/output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  double a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};
    
  const unsigned int dim=3;
  Tensor<2,dim> t(a);

  deallog << t.norm() << std::endl;
  t = 0;
  deallog << t.norm() << std::endl;

  Assert (t.norm() == 0, ExcInternalError());
  
  deallog << "OK" << std::endl;
}
