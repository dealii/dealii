//----------------------------  symmetric_tensor_02.cc  ---------------------------
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
//----------------------------  symmetric_tensor_02.cc  ---------------------------

// test symmetric 3x3 tensors

#include "../tests.h"
#include <base/symmetric_tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>

int main ()
{
  std::ofstream logfile("symmetric_tensor_02.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  SymmetricTensor<2,3> t;
  t[0][0] = 1;
  t[1][1] = 2;
  t[2][2] = 3;
  t[0][1] = 4;
  t[0][2] = 5;
  t[1][2] = 6;

  Assert (t[0][1] == t[1][0], ExcInternalError());

                                   // check that if a single element is
                                   // accessed, its transpose element gets the
                                   // same value
  t[1][0] = 14;
  Assert (t[0][1] == 14, ExcInternalError());

                                   // make sure transposition doesn't change
                                   // anything
  Assert (t == transpose(t), ExcInternalError());

                                   // check norm of tensor
  Assert (std::fabs(t.norm() - std::sqrt(1.*1+2*2+3*3+2*14*14+2*5*5+2*6*6))
          < 1e-14,
          ExcInternalError());

                                   // make sure norm is induced by scalar
                                   // product
  Assert (std::fabs (t.norm()*t.norm() - t*t) < 1e-14,
          ExcInternalError());

  deallog << "OK" << std::endl;
}
