//----------------------------  full_tensor_02.cc  ---------------------------
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
//----------------------------  full_tensor_02.cc  ---------------------------

// test full 3x3 tensors

#include "../tests.h"
#include <base/tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>

int main ()
{
  std::ofstream logfile("full_tensor_02.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  Tensor<2,3> t;
  t[0][0] = 1;
  t[1][1] = 2;
  t[2][2] = 3;
  t[0][1] = 4;
  t[1][0] = 4;
  t[0][2] = 5;
  t[2][0] = 5;
  t[1][2] = 6;
  t[2][1] = 6;

  Assert (t[0][1] == t[1][0], ExcInternalError());

                                   // make sure transposition doesn't change
                                   // anything
  Assert (t == transpose(t), ExcInternalError());

                                   // check norm of tensor
  Assert (std::fabs(t.norm() - std::sqrt(1.*1+2*2+3*3+2*4*4+2*5*5+2*6*6))
          < 1e-14,
          ExcInternalError());

  deallog << "OK" << std::endl;
}
