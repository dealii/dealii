// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// check PointerMatrix:checkTvmult

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

template <typename number>
void
checkTvmult(FullMatrix<number> &A,
            Vector<number> &    V,
            const std::string & name = "Test Matrix")
{
  deallog << "Tvmult" << std::endl;

  PointerMatrix<FullMatrix<number>, Vector<number>> P(&A, name.c_str());
  Vector<number>                                    O(A.m());
  P.Tvmult(O, V);

  // Check the dimensions of the result matrix
  DEAL_II_Assert(A.m() == O.size(), ExcInternalError());
  deallog << "Dimensions of result vector verified" << std::endl;

  // Verifying results with Method 2: O=A Transpose*V
  Vector<number> O_(A.m());
  A.Tvmult(O_, V);

  DEAL_II_Assert(O == O_, ExcInternalError());
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

  const double Adata[] = {2, 3, 4, 5};

  FullMatrix<double> A(2, 2);
  A.fill(Adata);

  Vector<double> V(2);
  V(0) = 1;
  V(1) = 2;

  checkTvmult<double>(A, V);
}
