//----------------------------  full_matrix_02.cc,v  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix_02.cc,v  ---------------------------

//check method mTmult of FullMatrix

#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>

#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/eigen.h>

const double entries_A[9] = { 1,2,3,4,5,6,7,8,9 };
const double entries_B[9] = { 2,1,1,1,2,3,2,1,2 };
const double entries_Z[9] = { 0,0,0,0,0,0,0,0,0 };

// Create a positive definite random matrix

int
main ()
{
  std::ofstream logfile("full_matrix_02/output");
  deallog << std::fixed;
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  std::srand(3391466);

  FullMatrix<double> A(3,3,entries_A);
  FullMatrix<double> B(3,3,entries_B);
  FullMatrix<double> Z(3,3,entries_Z);
  FullMatrix<double> C(3,3);
  FullMatrix<double> D(3,3);
 
  //compute C= A*B^T in two different ways and compare for equality
  Z.Tadd(1.,B);
  A.mmult(D,Z);
  A.mTmult(C,B);
  
  D.add(-1,C);
  Assert ( D.frobenius_norm() < 1e-15,
           ExcInternalError());
  deallog << "OK" << std::endl;
}    
