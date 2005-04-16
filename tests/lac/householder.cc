//--------------------------------------------------------------------
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
//--------------------------------------------------------------------

// Tests Householder class for QR-decomposition

#include "../tests.h"
#include <base/logstream.h>
#include <lac/householder.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>

#include <fstream>
#include <iostream>

const double rect[] =
{
      4., 3., 2., 1.,
      5., 8., 1., -2.,
      11., 13., -4., -5
};


int main()
{
  std::ofstream logfile("householder.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  FullMatrix<double> A(4,3,rect);
  Householder<double> H(A);
  
  Vector<double> u(4);
  Vector<double> v1(3);
  Vector<double> v2(3);  
  
  for (unsigned int i=0;i<u.size();++i)
    u(i) = i*i;
  deallog << "Distance " << H.least_squares(v1,u) << std::endl;
}
