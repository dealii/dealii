// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// Tests inverse SVD of LAPACKFullMatrix by comparing to vmult of FullMatrix

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

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


void test_rect(unsigned int m, unsigned int n, const double *values)
{
  std::ostringstream prefix;
  prefix << m << 'x' << n;
  deallog.push(prefix.str());

  FullMatrix<double> A(m,n,values);
  LAPACKFullMatrix<double> LA(m,n);
  LA = A;
  LA.compute_inverse_svd();

  deallog << "Singular values";
  for (unsigned int i=0; i<LA.n_rows(); ++i)
    deallog << ' ' << LA.singular_value(i);
  deallog << std::endl;

  Vector<double> u1(n);
  Vector<double> u2(n);
  Vector<double> v1(m);
  Vector<double> v2(m);

  for (unsigned int i=0; i<u1.size(); ++i)
    u1(i) = i*i;
  for (unsigned int i=0; i<v1.size(); ++i)
    v1(i) = i*i;

  // Test if LA is a left inverse of A
  A.vmult(v2,u1);
  LA.vmult(u2,v2);
  u2 -= u1;
  if (u2.l2_norm() < 1.e-12)
    deallog << "vmult ok" << std::endl;
  else
    deallog << "vmult error " << u2.l2_norm() << std::endl;

  // Test if LA is a right inverse of A
  LA.vmult(u2,v1);
  A.vmult(v2,u2);
  v2 -= v1;
  if (v2.l2_norm() < 1.e-12)
    deallog << "vmult ok" << std::endl;
  else
    deallog << "vmult error " << v2.l2_norm() << std::endl;

  // Test if LA^T is a left inverse
  // of A^T
  A.Tvmult(u2,v1);
  LA.Tvmult(v2, u2);
  v2 -= v1;
  if (v2.l2_norm() < 1.e-12)
    deallog << "Tvmult ok" << std::endl;
  else
    deallog << "Tvmult error " << v2.l2_norm() << std::endl;

  LA.Tvmult(v2,u1);
  A.Tvmult(u2,v2);
  u2 -= u1;
  if (u2.l2_norm() < 1.e-12)
    deallog << "Tvmult ok" << std::endl;
  else
    deallog << "Tvmult error " << u2.l2_norm() << std::endl;

  deallog.pop();
}

int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);

  test_rect(4,4,symm);
  test_rect(4,3,rect);
  test_rect(3,4,rect);
}
