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


//check method FullMatrix::triple_product

#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>

#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/eigen.h>

void diff(FullMatrix<double> &M)
{
  const double err = M.frobenius_norm();
  if (err < 1.e-14)
    deallog << "ok" << std::endl;
  else
    deallog << "oops " << err << std::endl;
}


void test (const unsigned int n, const unsigned int m)
{
  // Create some random matrices
  FullMatrix<double> A(n,n);
  FullMatrix<double> C(m,m);
  FullMatrix<double> B(n,m);
  FullMatrix<double> D(m,n);
  FullMatrix<double> Bt(m,n);
  FullMatrix<double> Dt(n,m);

  for (unsigned int i=0; i<n; ++i)
    {
      for (unsigned int j=0; j<n; ++j)
        A(i,j) = Testing::rand()/(1.0*RAND_MAX);
      for (unsigned int j=0; j<m; ++j)
        {
          B(i,j) = Bt(j,i) = Testing::rand()/(1.0*RAND_MAX);
          D(j,i) = Dt(i,j) = Testing::rand()/(1.0*RAND_MAX);
        }
    }
  for (unsigned int i=0; i<m; ++i)
    for (unsigned int j=0; j<m; ++j)
      C(i,j) = Testing::rand()/(1.0*RAND_MAX);

  // Compare first Schur complement
  // with mmult.
  FullMatrix<double> S1(m,m);
  S1.triple_product(A, D, B, false, false);

  FullMatrix<double> aux1(n,m);
  A.mmult(aux1,B);
  FullMatrix<double> aux2(m,m);
  D.mmult(aux2,aux1);
  aux2.add(-1., S1);
  diff(aux2);

  // Compare second Schur complement
  // with mmult
  FullMatrix<double> S2(n,n);
  S2.triple_product(C, B, D, false, false);

  FullMatrix<double> aux3(m,n);
  C.mmult(aux3,D);
  FullMatrix<double> aux4(n,n);
  B.mmult(aux4,aux3);
  aux4.add(-1., S2);
  diff(aux4);

  // Check transpose versions
  aux2 = 0.;
  aux2.triple_product(A, Dt, B, true, false);
  aux2.add(-1., S1);
  diff(aux2);
  aux2 = 0.;
  aux2.triple_product(A, D, Bt, false, true);
  aux2.add(-1., S1);
  diff(aux2);
  aux2 = 0.;
  aux2.triple_product(A, Dt, Bt, true, true);
  aux2.add(-1., S1);
  diff(aux2);

  aux4 = 0.;
  aux4.triple_product(C, Bt, D, true, false);
  aux4.add(-1., S2);
  diff(aux4);
  aux4 = 0.;
  aux4.triple_product(C, B, Dt, false, true);
  aux4.add(-1., S2);
  diff(aux4);
  aux4 = 0.;
  aux4.triple_product(C, Bt, Dt, true, true);
  aux4.add(-1., S2);
  diff(aux4);
}


int
main ()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);

  test(3,4);
  test(4,7);
}
