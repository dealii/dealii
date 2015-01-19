// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
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

// check PointerMatrix:checkConstructor2

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

template<typename number>
  void
  checkConstructor2(FullMatrix<number> &A)
  {
    deallog << "Init with matrix A" << std::endl;
    PointerMatrix<FullMatrix<number>, Vector<number> > P(&A);
    deallog << "Is matrix empty:" << P.empty() << std::endl;
  }

int
main()
{

  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const double Adata[] =
    { 2, 3, 4, 5 };

  FullMatrix<double> A(2, 2);
  A.fill(Adata);

  checkConstructor2<double>(A);
}
