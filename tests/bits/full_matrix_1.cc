// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// FullMatrix::copy_from could not be compiled if we copied from a
// sparse matrix. make sure this now works

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <fstream>


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const unsigned int N = 4;
  FullMatrix<double> f(N,N);

  SparsityPattern s (N,N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      s.add (i,j);
  s.compress ();

  SparseMatrix<double> sm(s);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      sm.set (i,j,i*j);

  f.copy_from (sm);

  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      {
        deallog << i << ' ' << j << ' ' << f(i,j) << std::endl;
        Assert (f(i,j) == sm(i,j), ExcInternalError());
      }

  deallog << "OK" << std::endl;
}
