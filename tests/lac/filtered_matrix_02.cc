// ---------------------------------------------------------------------
// $Id: filtered_matrix.cc 31349 2013-10-20 19:07:06Z maier $
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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
#include <deal.II/base/logstream.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

template<typename number>
void
checkTvmult(FullMatrix<number> &A, Vector<number> &V,
            bool expect_constrained_source = false)
{
  deallog << "Tvmult" << std::endl;

  FilteredMatrix < Vector<double> > F;
  F.initialize(A, expect_constrained_source);
  F.add_constraint(0, 1);

  Vector<number> O(A.n());

  F.Tvmult(O, V);

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
  { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  FullMatrix<double> A(3, 3);

  A.fill(Adata);

  Vector<double> V1(3);
  Vector<double> V2(3);

  V1(0) = V2(0) = 1;
  V1(1) = V2(1) = 2;
  V1(2) = V2(2) = 3;

  checkTvmult<double>(A, V1, false);
  checkTvmult<double>(A, V2, true);
}
