// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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


// Tests LAPACKFullMatrix::apply_lu_factorization in two different variants

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>



void test()
{
  const unsigned int n=11;
  LAPACKFullMatrix<double> A(n,n);
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<n; ++j)
      A(i,j) = (double)rand()/RAND_MAX;
  A.compute_lu_factorization();
  LAPACKFullMatrix<double> rhs_orig(n, 3);
  for (unsigned int i=0; i<3; ++i)
    for (unsigned int j=0; j<n; ++j)
      rhs_orig(j,i) = 2.32 * rand() / RAND_MAX - 0.9923;

  for (unsigned int transpose=0; transpose < 2; ++transpose)
    {
      LAPACKFullMatrix<double> rhs(rhs_orig);
      A.apply_lu_factorization(rhs, transpose);
      for (unsigned int i=0; i<3; ++i)
        {
          Vector<double> check(n);
          for (unsigned int j=0; j<n; ++j)
            check(j) = rhs_orig(j, i);
          A.apply_lu_factorization(check, transpose);
          for (unsigned int j=0; j<n; ++j)
            Assert(std::abs(check(j) - rhs(j,i)) < 1e-13, ExcInternalError());
        }
    }

  deallog << "OK" << std::endl;
}

int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test();
}
