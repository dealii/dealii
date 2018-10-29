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

// check PointerMatrix:checkConstructor2

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

template <typename number>
void
checkConstructor2(FullMatrix<number> &A)
{
  deallog << "Init with matrix A" << std::endl;
  PointerMatrix<FullMatrix<number>, Vector<number>> P(&A);
  deallog << "Is matrix empty:" << P.empty() << std::endl;
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

  checkConstructor2<double>(A);
}
