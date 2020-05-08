// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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


// Similar to precondition_chebyshev_01 but using a separate preconditioner
// class for providing the inverse diagonal that goes through another code
// path in PreconditionChebyshev


#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"



class DiagonalMatrixManual
{
public:
  DiagonalMatrixManual()
  {}

  void
  set_vector_one(const unsigned int size)
  {
    diagonal.reinit(size);
    diagonal = 1;
  }

  void
  reinit(const FullMatrix<double> f)
  {
    diagonal.reinit(f.m());
    for (unsigned int i = 0; i < f.m(); ++i)
      diagonal(i) = 1. / f(i, i);
    deallog.get_file_stream() << "DEAL::";
    diagonal.print(deallog.get_file_stream());
  }

  typename Vector<double>::size_type
  m() const
  {
    return diagonal.size();
  }

  void
  vmult(Vector<double> &dst, const Vector<double> &src) const
  {
    dst = src;
    dst.scale(diagonal);
  }

private:
  Vector<double> diagonal;
};


void
check()
{
  const unsigned int size = 10;
  FullMatrix<double> m(size, size);
  for (unsigned int i = 0; i < size; ++i)
    m(i, i) = i + 1;

  Vector<double> in(size), out(size);
  for (unsigned int i = 0; i < size; ++i)
    in(i) = random_value<double>();

  PreconditionChebyshev<FullMatrix<double>,
                        Vector<double>,
                        DiagonalMatrixManual>
                                                              prec;
  PreconditionChebyshev<FullMatrix<double>,
                        Vector<double>,
                        DiagonalMatrixManual>::AdditionalData data;
  data.smoothing_range = 2 * size;
  data.degree          = 4;
  data.preconditioner.reset(new DiagonalMatrixManual());
  data.preconditioner->set_vector_one(size);
  prec.initialize(m, data);

  deallog << "Exact inverse:     ";
  for (unsigned int i = 0; i < size; ++i)
    deallog << in(i) / m(i, i) << " ";
  deallog << std::endl;

  deallog << "Check  vmult ones: ";
  prec.vmult(out, in);
  for (unsigned int i = 0; i < size; ++i)
    deallog << out(i) << " ";
  deallog << std::endl;

  deallog << "Check Tvmult ones: ";
  prec.Tvmult(out, in);
  for (unsigned int i = 0; i < size; ++i)
    deallog << out(i) << " ";
  deallog << std::endl;

  data.preconditioner->reinit(m);
  data.smoothing_range = 2;
  data.degree          = 5;
  prec.initialize(m, data);

  deallog << "Check  vmult diag: ";
  prec.vmult(out, in);
  for (unsigned int i = 0; i < size; ++i)
    deallog << out(i) << " ";
  deallog << std::endl;

  deallog << "Check Tvmult diag: ";
  prec.Tvmult(out, in);
  for (unsigned int i = 0; i < size; ++i)
    deallog << out(i) << " ";
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
