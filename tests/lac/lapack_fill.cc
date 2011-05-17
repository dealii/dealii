//-------------------------  lapack_fill.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------  lapack_fill.cc  ---------------------------


#include "../tests.h"
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <fstream>

// A.fill() produced an ExcIndexRange(r,0,m()) exception with
// the additional Information: Index 6 is not in [0,3[.
// Bug reported by Florian Prill

int main () 
{
  std::ofstream logfile("lapack_fill/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
				   // matrix sizes
  const unsigned int m = 3;
  const unsigned int n = 10;       

  LAPACKFullMatrix<double> A(n);
  FullMatrix<double>       C(m);
				   // fill some entries:
  C(0,0)     = 1.0;
  C(m-1,m-1) = 1.0;
				   // insert C into A's middle:
  A.fill(C,
	 3,3,
	 0,0);
				   // check some values
  Assert(A(3,3)==1, ExcInternalError());
  Assert(A(5,5)==1, ExcInternalError());
  
  deallog << "OK" << std::endl;
}
