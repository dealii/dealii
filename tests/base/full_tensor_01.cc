//----------------------------  full_tensor_01.cc  ---------------------------
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
//----------------------------  full_tensor_01.cc  ---------------------------

// test full 2x2 tensors

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("full_tensor_01/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  Tensor<2,2> t;
  t[0][0] = 1;
  t[1][1] = 2;
  t[0][1] = 4;
  t[1][0] = 4;
                                   // make sure transposition doesn't change
                                   // anything
  Assert (t == transpose(t), ExcInternalError());

                                   // check norm of tensor
  Assert (std::fabs(t.norm() - std::sqrt(1.*1+2*2+2*4*4)) < 1e-14,
          ExcInternalError());

  deallog << "OK" << std::endl;
}
