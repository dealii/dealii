//----------------------------  symmetric_tensor_01.cc  ---------------------------
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
//----------------------------  symmetric_tensor_01.cc  ---------------------------

// test symmetric 2x2 tensors

#include "../tests.h"
#include <base/symmetric_tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>

int main ()
{
  std::ofstream logfile("symmetric_tensor_01.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  SymmetricTensor<2,2> t;
  t[0][0] = 1;
  t[1][1] = 2;
  t[0][1] = 3;

  Assert (t[0][1] == t[1][0], ExcInternalError());

                                   // check that if a single element is
                                   // accessed, its transpose element gets the
                                   // same value
  t[1][0] = 4;
  Assert (t[0][1] == 4, ExcInternalError());

                                   // make sure transposition doesn't change
                                   // anything
  Assert (t == transpose(t), ExcInternalError());

                                   // check norm of tensor
  Assert (std::fabs(t.norm() - std::sqrt(1.*1+2*2+2*4*4)) < 1e-14,
          ExcInternalError());

                                   // make sure norm is induced by scalar
                                   // product
  Assert (std::fabs (t.norm()*t.norm() - t*t) < 1e-14,
          ExcInternalError());

  deallog << "OK" << std::endl;
}
