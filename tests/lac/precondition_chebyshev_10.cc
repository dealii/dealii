// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Tests PreconditionChebyshev::vmult_with_last_residual_norm() and
// PreconditionChebyshev::step_with_last_residual_norm()


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
check(const unsigned int size)
{
  deallog << "Check operation of size " << size << std::endl;
  FullMatrixModified m(size, size);
  for (unsigned int i = 0; i < size; ++i)
    m(i, i) = i + 1;

  Vector<double> in(size), out(size), tmp(size);
  for (unsigned int i = 0; i < size; ++i)
    in(i) = random_value<double>();

  PreconditionChebyshev<FullMatrixModified, Vector<double>> prec;
  PreconditionChebyshev<FullMatrixModified, Vector<double>>::AdditionalData
    data;
  data.smoothing_range      = 1.5 * size;
  data.eig_cg_n_iterations  = 0;
  data.max_eigenvalue       = size + 1;
  data.degree               = 8;
  data.eigenvalue_algorithm = internal::EigenvalueAlgorithm::power_iteration;
  prec.initialize(m, data);

  const double res_norm = prec.vmult_with_last_residual_norm(out, in);
  deallog << "Computed residual norm vmult(): " << res_norm << std::endl;

  prec.set_degree(data.degree - 1);
  prec.vmult(out, in);

  m.vmult(tmp, out);
  tmp.sadd(-1.0, 1.0, in);
  deallog << "Actual residual norm vmult():   " << tmp.l2_norm() << std::endl;

  const double res_norm_step = prec.step_with_last_residual_norm(out, in);
  deallog << "Computed residual norm step():  " << res_norm_step << std::endl;

  prec.set_degree(data.degree - 1);
  prec.vmult(out, in);
  prec.set_degree(data.degree - 2);
  prec.step(out, in);

  m.vmult(tmp, out);
  tmp.sadd(-1.0, 1.0, in);
  deallog << "Actual residual norm step():    " << tmp.l2_norm() << std::endl;
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(10);
  deallog.attach(logfile);

  MultithreadInfo::set_thread_limit(1);

  check(10);
  check(16);
  // Check a few sizes to ensure vectorization granularity is kept correct
  for (unsigned int i = 23; i < 34; ++i)
    check(i);

  return 0;
}
