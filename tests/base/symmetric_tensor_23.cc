//----------------------------  symmetric_tensor_23.cc  ---------------------------
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
//----------------------------  symmetric_tensor_23.cc  ---------------------------

// check operator<< for SymmetricTensor<2,dim> and SymmetricTensor<4,dim>

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("symmetric_tensor_23/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  {
    SymmetricTensor<2,1> t;
    t[0][0] = 1;

    double x[1] = { 1 };
    Assert ((t == SymmetricTensor<2,1>(x)),
	    ExcInternalError());
  }

  {
    SymmetricTensor<2,2> t;
    t[0][0] = 1;
    t[1][1] = 2;
    t[0][1] = 3;

    double x[3] = { 1, 2, 3 };
    Assert ((t == SymmetricTensor<2,2>(x)),
	    ExcInternalError());
  }

  {
    SymmetricTensor<2,3> t;
    t[0][0] = 1;
    t[1][1] = 2;
    t[2][2] = 3;
    t[0][1] = 4;
    t[0][2] = 5;
    t[1][2] = 6;

    double x[6] = { 1, 2, 3, 4, 5, 6 };
    Assert ((t == SymmetricTensor<2,3>(x)),
	    ExcInternalError());
  }
  
  deallog << "OK" << std::endl;
}
