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

// Check ProductMatrix::vmult_add

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>
#include <cmath>

template<typename number>
  void
  checkVmult_add(FullMatrix<number> &A, FullMatrix<number> &B,
      Vector<number> &V)
  {
    deallog << "vmult_add" << std::endl;
    GrowingVectorMemory < Vector<number> > mem;
    ProductMatrix < Vector<number> > P(A, B, mem);

    deallog
        << "Result vector set to all ones and to be added with product result vector"
        << std::endl;
    Vector<number> O(V.size());
    for (unsigned int i = 0; i < O.size(); ++i)
      O(i) = 1;

    P.vmult_add(O, V);

    // Check the dimensions of the result matrix
    Assert(A.m() == O.size(), ExcInternalError());
    deallog << "Dimensions of result vector verified" << std::endl;

    // Verifying results with Method 2: O=(A*B)*V
    FullMatrix<number> AB_Product(A.m(), B.n());
    Vector<number> O_(A.m());
    for (unsigned int i = 0; i < O_.size(); ++i)
      O_(i) = 1;

    A.mmult(AB_Product, B, false);
    AB_Product.vmult_add(O_, V);

    Assert(O == O_, ExcInternalError());
    deallog << "Result vector data verified" << std::endl;

    for (unsigned int i = 0; i < O.size(); ++i)
      deallog << O(i) << '\t';
    deallog << std::endl;
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

  const double Bdata[] =
    { 0, 1, 2, 3 };

  FullMatrix<double> A(2, 2);
  FullMatrix<double> B(2, 2);

  A.fill(Adata);
  B.fill(Bdata);

  Vector<double> V(2);
  V(0) = 1;
  V(1) = 2;

  checkVmult_add<double>(A, B, V);

}
