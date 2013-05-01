//----------------------------  full_matrix_07.cc,v  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix_07.cc,v  ---------------------------

//check method mTmult of FullMatrix, symmetric case

#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>

#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>

const double entries_A[9] = { 1,2,3,4,5,6,7,8,9 };
const double compare[9] = { 14,32,50,32,77,122,50,122,194 };

int
main ()
{
  std::ofstream logfile("full_matrix_07/output");
  deallog << std::fixed;
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  FullMatrix<double> A(3,3,entries_A);
  FullMatrix<double> C(3,3);
  FullMatrix<double> D(3,3,compare);
 
  //compute C= A*A^T
  A.mTmult(C,A);

  C.add(-1., D);
  Assert(C.frobenius_norm() < 1e-12, ExcInternalError());
  
  deallog << "OK" << std::endl;
}    
