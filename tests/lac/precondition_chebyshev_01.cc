// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests PreconditionChebyshev::vmult and PreconditionChebyshev::Tvmult


#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"



class FullMatrixModified : public FullMatrix<double>
{
public:
  FullMatrixModified(unsigned int size1, unsigned int size2)
    : FullMatrix<double>(size1, size2)
  {}

  double
  el(unsigned int i, unsigned int j) const
  {
    return this->operator()(i, j);
  }
};


void
check()
{
  const unsigned int size = 10;
  FullMatrixModified m(size, size);
  for (unsigned int i = 0; i < size; ++i)
    m(i, i) = i + 1;

  Vector<double> in(size), out(size);
  for (unsigned int i = 0; i < size; ++i)
    in(i) = random_value<double>();

  PreconditionChebyshev<FullMatrixModified, Vector<double>> prec;
  PreconditionChebyshev<FullMatrixModified, Vector<double>>::AdditionalData
    data;
  data.smoothing_range = 2 * size;
  data.degree          = 4;
  prec.initialize(m, data);

  deallog << "Exact inverse:     ";
  for (unsigned int i = 0; i < size; ++i)
    deallog << in(i) / m(i, i) << ' ';
  deallog << std::endl;

  deallog << "Check  vmult orig: ";
  prec.vmult(out, in);
  for (unsigned int i = 0; i < size; ++i)
    deallog << out(i) << ' ';
  deallog << std::endl;

  deallog << "Check Tvmult orig: ";
  prec.Tvmult(out, in);
  for (unsigned int i = 0; i < size; ++i)
    deallog << out(i) << ' ';
  deallog << std::endl;

  Vector<double> matrix_diagonal(size);
  matrix_diagonal     = 1;
  auto preconditioner = std::make_shared<DiagonalMatrix<Vector<double>>>();
  preconditioner->reinit(matrix_diagonal);
  data.preconditioner = std::move(preconditioner);
  prec.initialize(m, data);

  deallog << "Check  vmult diag: ";
  prec.vmult(out, in);
  for (unsigned int i = 0; i < size; ++i)
    deallog << out(i) << ' ';
  deallog << std::endl;

  deallog << "Check Tvmult diag: ";
  prec.Tvmult(out, in);
  for (unsigned int i = 0; i < size; ++i)
    deallog << out(i) << ' ';
  deallog << std::endl;
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
