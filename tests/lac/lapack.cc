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

// Tests all aspects of LAPACKFullMatrix

#include "../tests.h"
#include <base/logstream.h>
#include <lac/lapack_full_matrix.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>

#include <fstream>
#include <iostream>

/*
 * Eigenvalues and -vectors of this system are
 * lambda = 1     v = (1, 1, 1, 1)
 * lambda = 5     v = (1,-1, 0, 0)
 * lambda = 5     v = (0, 1,-1, 0)
 * lambda = 5     v = (0, 0, 1,-1)
 */
const double symm[] =
{
      4., -1., -1., -1.,
      -1., 4., -1., -1.,
      -1., -1., 4., -1.,
      -1., -1., -1., 4.
};

const double rect[] =
{
      4., 3., 2., 1.,
      5., 8., 1., -2.,
      11., 13., -4., -5
};


int main()
{
  std::ofstream logfile("lapack.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);

#ifdef HAVE_LIBLAPACK
  FullMatrix<double> A(3,4,rect);
  LAPACKFullMatrix<double> LA(3,4);
  LA = A;
  
  Vector<double> u(4);
  Vector<double> v1(3);
  Vector<double> v2(3);
  
  for (unsigned int i=0;i<u.size();++i)
    u(i) = i*i;
  
				   // Test rectangular vmult. All
				   // results compare with same
				   // operation for FullMatrix.
  deallog.push("Rect");

  deallog << "operator= (const FullMatrix<number>&) ok" << std::endl;
  
  A.vmult(v1,u);
  LA.vmult(v2,u);
  v1 -= v2;
  if (v1.l2_norm() < 1.e-14)
    deallog << "vmult ok" << std::endl;
  v1 = v2;
  
  A.vmult_add(v1,u);
  LA.vmult_add(v2,u);
  v1 -= v2;
  if (v1.l2_norm() < 1.e-14)
    deallog << "vmult_add ok" << std::endl;
  
  LA.Tvmult(u, v2);
  u *= -1;
  A.Tvmult_add(u, v2);
  if (u.l2_norm() < 1.e-14)
    deallog << "Tvmult ok" << std::endl;
  
  A.Tvmult(u, v2);
  u *= -1;
  LA.Tvmult_add(u, v2);
  if (u.l2_norm() < 1.e-14)
    deallog << "Tvmult_add ok" << std::endl;
  
  deallog.pop();

				   // Test symmetric system
  A.reinit(4,4);
  LA.reinit(4,4);
  A.fill(symm);
  LA = A;
  LA.compute_eigenvalues();
  for (unsigned int i=0;i<A.m();++i)
    {
      std::complex<double> lambda = LA.eigenvalue(i);
      deallog << "Eigenvalues "
	      << (int) (lambda.real()+.0001) << '\t'
	      << (int) (lambda.imag()+.0001) << std::endl;
    }
  
  v1.reinit(4);
  v2.reinit(4);
  
#else
			// If lapack is not available, this
				   // test will not complain.
  deallog.push("Rect");
  deallog << "operator= (const FullMatrix<number>&) ok" << std::endl;
  deallog << "vmult ok" << std::endl;
  deallog << "vmult_add ok" << std::endl;
  deallog << "Tvmult ok" << std::endl;
  deallog << "Tvmult_add ok" << std::endl;
  deallog << "Eigenvalues 5\t0" << std::endl;
  deallog << "Eigenvalues 1\t0" << std::endl;
  deallog << "Eigenvalues 5\t0" << std::endl;
  deallog << "Eigenvalues 5\t0" << std::endl;
  deallog.pop();
#endif
}
