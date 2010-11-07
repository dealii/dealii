//--------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------

// LAPACKFullMatrix::compute_lu_factorization
// LAPACKFullMatrix::apply_lu_factorization

#include "../tests.h"
#include <base/logstream.h>
#include <lac/lapack_full_matrix.h>
#include <lac/vector.h>

#include <fstream>
#include <iostream>


// Fill a matrix with the values of the Hilbert matrix
template <typename number>
void
hilbert (LAPACKFullMatrix<number>& M,
	 const bool nonsymmetric)
{
  const unsigned int n = M.n_rows();
  for (unsigned int i=0;i<n;++i)
    for (unsigned int j=0;j<n;++j)
      M(i,j) = 1./(i+j+1.)
	       * (nonsymmetric
		  ? ((i<j) ? -1. : 1.)
		  : 1.);
}


// Multiply some vectors with the matrix and its transpose, then
// compute and apply LU factorization and see if the results are equal
// to the original vector.

void
test(const unsigned int size, const bool nonsymmetric)
{
  LAPACKFullMatrix<double> M(size, size);
  hilbert(M, nonsymmetric);
  
  Vector<double> u(size);
  Vector<double> v(size);
  Vector<double> x(size);
  Vector<double> y(size);
  
  for (unsigned int i=0;i<size;++i)
    {
      u(i) = i+2.;
      x(i) = i+2.;
    }
  M.vmult(v,u);
  M.Tvmult(y,x);
  M.compute_lu_factorization();
  M.apply_lu_factorization(v, false);
  M.apply_lu_factorization(y, true);

  v -= u;
  y -= x;

  deallog << v.l2_norm() << std::endl;
  deallog << y.l2_norm() << std::endl;
}


int main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test(4, false);
  test(4, true);
  test(5, false);
  test(5, true);
}
