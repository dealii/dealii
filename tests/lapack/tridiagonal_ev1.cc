// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


// Tests compute_eigenvalues() and eigenvalue() of TridiagonalMatrix

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/tridiagonal_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>


// Build the one dimensional discrete Laplacian for n intervals
template <typename number>
void test_laplacian(unsigned int n)
{
  TridiagonalMatrix<number> M(n-1, true);
  for (unsigned int i=0; i<n-2; ++i)
    {
      M(i,i) = 2.;
      M(i,i+1) = -1.;
    }
  M(n-2,n-2) = 2.;

  M.compute_eigenvalues();
  for (unsigned int i=0; i<5; ++i)
    deallog << '\t' << M.eigenvalue(i)*n *n;
  deallog << "\t cond " << M.eigenvalue(n-2)/M.eigenvalue(0) << std::endl;
}


int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_laplacian<double>(10);
  test_laplacian<double>(20);
  test_laplacian<double>(40);
  test_laplacian<double>(80);
}
