// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests DiagonalMatrix


#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"



void
check()
{
  const unsigned int size = 10;
  Vector<double>     vec(size);
  for (unsigned int i = 0; i < size; ++i)
    vec(i) = 2;
  DiagonalMatrix<Vector<double>> mat;
  mat.reinit(vec);

  Vector<double> in(size), out(size), exact(size);
  for (unsigned int i = 0; i < size; ++i)
    in(i) = 0.5;
  exact = 1;

  mat.vmult(out, in);
  out -= exact;
  deallog << "Error vmult: " << out.linfty_norm() << std::endl;

  mat.Tvmult(out, in);
  out -= exact;
  deallog << "Error Tvmult: " << out.linfty_norm() << std::endl;

  out = exact;
  mat.vmult_add(out, in);
  out -= exact;
  out -= exact;
  deallog << "Error vmult_add: " << out.linfty_norm() << std::endl;

  out = exact;
  mat.Tvmult_add(out, in);
  out -= exact;
  out -= exact;
  deallog << "Error Tvmult_add: " << out.linfty_norm() << std::endl;

  for (unsigned int i = 0; i < size; ++i)
    {
      mat(i, i) = i + 1;
      in(i)     = 1. / mat(i, i);
    }

  mat.vmult(out, in);
  out -= exact;
  deallog << "Error vmult set 1: " << out.linfty_norm() << std::endl;

  for (unsigned int i = 0; i < size; ++i)
    {
      mat.get_vector()(i) = 1. / (i + 1);
      in(i)               = i + 1;
    }
  mat.vmult(out, in);
  out -= exact;
  deallog << "Error vmult set 2: " << out.linfty_norm() << std::endl;

  for (unsigned int i = 0; i < size; ++i)
    {
      const double mat_entry =
        const_cast<const DiagonalMatrix<Vector<double>> &>(mat)(i, i);
      mat(i, i) = in(i);
      in(i)     = mat_entry;
    }
  mat.vmult(out, in);
  out -= exact;
  deallog << "Error vmult set 3: " << out.linfty_norm() << std::endl;
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);

  check();

  return 0;
}
