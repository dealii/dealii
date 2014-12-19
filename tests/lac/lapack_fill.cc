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



#include "../tests.h"
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <fstream>

// A.fill() produced an ExcIndexRange(r,0,m()) exception with
// the additional Information: Index 6 is not in [0,3[.
// Bug reported by Florian Prill

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  // matrix sizes
  const unsigned int m = 3;
  const unsigned int n = 10;

  LAPACKFullMatrix<double> A(n);
  FullMatrix<double>       C(m);
  // fill some entries:
  C(0,0)     = 1.0;
  C(m-1,m-1) = 1.0;
  // insert C into A's middle:
  A.fill(C,
         3,3,
         0,0);
  // check some values
  Assert(A(3,3)==1, ExcInternalError());
  Assert(A(5,5)==1, ExcInternalError());

  deallog << "OK" << std::endl;
}
