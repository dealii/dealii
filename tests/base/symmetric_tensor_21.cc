//----------------------------  symmetric_tensor_21.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  symmetric_tensor_21.cc  ---------------------------

// check SymmetricTensor<4,dim>::operator= (double)

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("symmetric_tensor_21/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  const unsigned int dim=3;
  SymmetricTensor<4,dim> t;
  t[0][0][0][0] = t[1][0][1][0] = t[1][1][1][1]
                = t[2][2][2][2] = t[2][0][2][0] = 3;

  deallog << t.norm() << std::endl;
  t = 0;
  deallog << t.norm() << std::endl;

  Assert (t.norm() == 0, ExcInternalError());
  
  deallog << "OK" << std::endl;
}
