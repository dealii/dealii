//---------------------------------------------------------------------------
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
//---------------------------------------------------------------------------

// test symmetric 2x2x2x2 tensors

#include <base/symmetric_tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>

int main ()
{
  std::ofstream logfile("symmetric_tensor_03.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  SymmetricTensor<4,2> t;
  t[0][0][0][0] = 1;
  t[1][1][1][1] = 2;
  t[0][1][0][1] = 3;

  Assert (t[0][1][0][1] == t[1][0][1][0], ExcInternalError());

                                   // check that if a single element is
                                   // accessed, its transpose element gets the
                                   // same value
  t[1][0][0][1] = 4;
  Assert (t[0][1][1][0] == 4, ExcInternalError());

                                   // make sure transposition doesn't change
                                   // anything
  Assert (t == transpose(t), ExcInternalError());

                                   // check norm of tensor
  deallog << t.norm() << std::endl;

                                   // make sure norm is induced by scalar
                                   // product
  double norm_sqr = 0;
  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<2; ++j)
      for (unsigned int k=0; k<2; ++k)
	for (unsigned int l=0; l<2; ++l)
	  norm_sqr += t[i][j][k][l] * t[i][j][k][l];
  
  Assert (std::fabs (t.norm()*t.norm() - norm_sqr) < 1e-14,
          ExcInternalError());

  deallog << "OK" << std::endl;
}
