// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// check SparseMatrix::mmult. this function has a default argument
// that could previously not be instantiated because it was in a
// non-deduced context but that should not be possible to omit
//
// this test checks SparseMatrix::mmult with additional arguments

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
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<n; ++j)
      {
        sp.add (i,j);
        C_sp.add (i,j);
      }
  sp.compress ();
  C_sp.compress ();

  SparseMatrix<double> A(sp);
  SparseMatrix<double> B(sp);
  SparseMatrix<double> C(C_sp);

  for (unsigned int i=0; i<n; ++i)
    {
      for (unsigned int j=0; j<n; ++j)
        A.set(i,j,Testing::rand());
      for (unsigned int j=0; j<n; ++j)
        B.set(i,j,Testing::rand());
    }

  Vector<double> v(n);
  for (unsigned int j=0; j<n; ++j)
    v(j) = Testing::rand();

  // now form the matrix-matrix product and
  // initialize a test rhs
  A.mmult (C, B, v);

  Vector<double> x(n), y(n), z(n), tmp(n);
  for (unsigned int j=0; j<n; ++j)
    x(j) = Testing::rand();

  // then test for correctness
  C.vmult (y, x);

  B.vmult (tmp, x);
  tmp.scale (v);
  A.vmult (z, tmp);

  y -= z;
  Assert (y.l2_norm() <= 1e-12 * z.l2_norm(),
          ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main ()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  Testing::srand(3391466);

  test(3);
  test(7);
}
