//--------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2010 by the deal.II authors
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


// Build the one dimensional discrete Laplacian for n intervals
template <typename number>
void test_laplacian(unsigned int n)
{
  TridiagonalMatrix<number> M(n-1, true);
  for (unsigned int i=0;i<n-2;++i)
    {
      M(i,i) = 2.;
      M(i,i+1) = -1.;
    }
  M(n-2,n-2) = 2.;

  M.compute_eigenvalues();
  for (unsigned int i=0;i<5;++i)
    deallog << '\t' << M.eigenvalue(i)*n*n;
  deallog << "\t cond " << M.eigenvalue(n-2)/M.eigenvalue(0) << std::endl;
}


int main()
{ 
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_laplacian<double>(10);
  test_laplacian<double>(20);
  test_laplacian<double>(40);
  test_laplacian<double>(80);
}
