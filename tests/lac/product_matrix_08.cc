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

// Check ProductMatrix::clear

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
  checkClear(FullMatrix<number> &A, FullMatrix<number> &B,
      FullMatrix<number> &C, FullMatrix<number> &D)
  {
    deallog << "clear" << std::endl;
    deallog << "Init with vector memory, matrix 1 and matrix 2" << std::endl;
    GrowingVectorMemory < Vector<number> > mem;
    ProductMatrix < Vector<number> > P(A, B, mem);

    deallog << "Multiplying with all ones vector" << std::endl;
    Vector<number> V(B.n());
    for (unsigned int i = 0; i < V.size(); ++i)
      V(i) = 1;

    Vector<number> O(A.m());

    P.vmult(O, V);

    // Check the dimensions of the result vector
    Assert(A.m() == O.size(), ExcInternalError());
    deallog << "Dimensions of result vector verified" << std::endl;

    // Verifying results with Method 2: O=(A*B)*V
    FullMatrix<number> AB_Product(A.m(), B.n());
    Vector<number> O_(A.m());

    A.mmult(AB_Product, B, false);
    AB_Product.vmult(O_, V);

    Assert(O == O_, ExcInternalError());
    deallog << "Result vector data verified" << std::endl;

    for (unsigned int i = 0; i < O.size(); ++i)
      deallog << O(i) << '\t';
    deallog << std::endl;

    deallog << "Clearing product matrix" << std::endl;
    P.clear();

    deallog << "Reinitializing product matrix with matrix 3 and matrix 4"
        << std::endl;
    P.reinit(C, D);

    deallog << "Multiplying with all ones vector" << std::endl;
    Vector<number> _V(D.n());
    for (unsigned int i = 0; i < _V.size(); ++i)
      _V(i) = 1;

    Vector<number> _O(C.m());

    P.vmult(_O, _V);

    // Check the dimensions of the result vector
    Assert(C.m() == _O.size(), ExcInternalError());
    deallog << "Dimensions of result vector verified" << std::endl;

    // Verifying results with Method 2: _O=(C*D)*_V
    FullMatrix<number> CD_Product(C.m(), D.n());
    Vector<number> _O_(C.m());

    C.mmult(CD_Product, D, false);
    CD_Product.vmult(_O_, _V);

    Assert(_O == _O_, ExcInternalError());
    deallog << "Result vector data verified" << std::endl;

    for (unsigned int i = 0; i < _O.size(); ++i)
      deallog << _O(i) << '\t';
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

  const double Cdata[] =
    { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  const double Ddata[] =
    { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

  FullMatrix<double> A(2, 2);
  FullMatrix<double> B(2, 2);
  FullMatrix<double> C(3, 3);
  FullMatrix<double> D(3, 3);

  A.fill(Adata);
  B.fill(Bdata);
  C.fill(Cdata);
  D.fill(Ddata);

  checkClear<double>(A, B, C, D);

}
