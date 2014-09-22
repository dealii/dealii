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


// Test vmult and Tvmult of PointerMatrix and TransposeMatrix

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/transpose_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  FullMatrix<double> A(5,5);
  unsigned int k=0;
  for (unsigned int i=0; i<A.m(); ++i)
    for (unsigned int j=0; j<A.n(); ++j)
      A(i,j) = ++k;

  PointerMatrix<FullMatrix<double>, Vector<double> > P(&A, "P");
  TransposeMatrix<FullMatrix<double>, Vector<double> > T(&A, "T");

  Vector<double> x(5);
  Vector<double> y(5);
  Vector<double> y2(5);
  Vector<double> diff(5);

  for (unsigned int i=0; i<x.size(); ++i)
    {
      x = 0.;
      x(i) = 1.;
      A.vmult(y,x);
      P.vmult(y2,x);
      diff = y;
      diff.add(-1., y2);
      deallog << "P.vmult:  diff " << diff.l2_norm() << std::endl;
      T.Tvmult(y2,x);
      diff = y;
      diff.add(-1., y2);
      deallog << "T.Tvmult: diff " << diff.l2_norm() << std::endl;

      A.Tvmult(y,x);
      P.Tvmult(y2,x);
      diff = y;
      diff.add(-1., y2);
      deallog << "P.Tvmult: diff " << diff.l2_norm() << std::endl;
      T.vmult(y2,x);
      diff = y;
      diff.add(-1., y2);
      deallog << "T.vmult:  diff " << diff.l2_norm() << std::endl;


    }
}
