//--------------------------------------------------------------------
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
//--------------------------------------------------------------------

// Tests compute_eigenvalues() and eigenvalue() of TridiagonalMatrix

#include "../tests.h"
#include <base/logstream.h>
#include <lac/tridiagonal_matrix.h>
#include <lac/vector.h>

#include <fstream>
#include <iostream>


template <typename number>
void test_laplacian(unsigned int n)
{
  TridiagonalMatrix<number> M(n, true);
  for (unsigned int i=0;i<n-1;++i)
    {
      M(i,i) = 2.;
      M(i,i+1) = -1.;
    }
  M(n-1,n-1) = 2.;

  M.compute_eigenvalues();
  for (unsigned int i=0;i<n;++i)
    deallog << i << '\t' << M.eigenvalue(i);
}


int main()
{
  std::ofstream logfile("tridiagonal_ev1/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_laplacian<double>(50);
}
