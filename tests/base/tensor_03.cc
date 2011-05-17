//----------------------------  tensor_03.cc  ---------------------------
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
//----------------------------  tensor_03.cc  ---------------------------

// check Tensor<1,dim>::operator= (double)

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("tensor_03/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  double a[3] = {1, 2, 3};
    
  const unsigned int dim=3;
  Tensor<1,dim> t(a);

  deallog << t.norm() << std::endl;
  t = 0;
  deallog << t.norm() << std::endl;

  Assert (t.norm() == 0, ExcInternalError());
  
  deallog << "OK" << std::endl;
}
