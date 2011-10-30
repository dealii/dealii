//----------------------------  sparse_matrix_mmult_01.cc,v  ---------------------------
//    $Id$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_matrix_mmult_01.cc,v  ---------------------------

// check SparseMatrix::mmult. this function has a default argument
// that could previously not be instantiated because it was in a
// non-deduced context but that should not be possible to omit
//
// this test checks SparseMatrix::mmult without additional arguments

#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>

#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>


void test (const unsigned int n)
{
				   // Create some random full matrices in the
				   // data structures of a sparse matrix
  SparsityPattern sp (n,n);
  SparsityPattern C_sp (n,n);
  for (unsigned int i=0;i<n;++i)
    for (unsigned int j=0;j<n;++j)
      {
	sp.add (i,j);
	C_sp.add (i,j);
      }
  sp.compress ();
  C_sp.compress ();

  SparseMatrix<double> A(sp);
  SparseMatrix<double> B(sp);
  SparseMatrix<double> C(C_sp);

  for (unsigned int i=0;i<n;++i)
    {
      for (unsigned int j=0;j<n;++j)
	A.set(i,j,std::rand());
      for (unsigned int j=0;j<n;++j)
	B.set(i,j,std::rand());
    }

				   // now form the matrix-matrix product and
				   // initialize a test rhs
  A.mmult (C, B);

  Vector<double> x(n), y(n), z(n), tmp(n);
  for (unsigned int j=0;j<n;++j)
    x(j) = std::rand();

				   // then test for correctness
  C.vmult (y, x);

  B.vmult (tmp, x);
  A.vmult (z, tmp);

  y -= z;
  Assert (y.l2_norm() <= 1e-12 * z.l2_norm(),
	  ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main ()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  std::srand(3391466);

  test(3);
  test(7);
}
