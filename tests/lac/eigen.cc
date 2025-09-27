// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/lac/eigen.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"

#include "../testmatrix.h"

int
main()
{
  std::ofstream logfile("output");
  //  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);

  GrowingVectorMemory<> mem;
  SolverControl         control(1000, 1.e-5);

  const unsigned int size = 10;
  const unsigned int dim  = (size - 1) * (size - 1);

  /*
   * Compute minimal and maximal
   * eigenvalues of the 5-point
   * stencil matrix
   * (Hackbusch:Iterative Loesung...,
   * Satz 4.1.1)
   */
  const double h          = 1. / size;
  const double s          = std::sin(numbers::PI * h / 2.);
  const double c          = std::cos(numbers::PI * h / 2.);
  const double lambda_max = 8. * c * c;
  const double lambda_min = 8. * s * s;

  FDMatrix        testproblem(size, size);
  SparsityPattern structure(dim, dim, 5);
  testproblem.five_point_structure(structure);
  structure.compress();
  SparseMatrix<double> A(structure);
  testproblem.five_point(A);

  Vector<double> u(dim);
  u = 1.;
  double lambda;

  EigenPower<> mises(control, mem);
  mises.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << " Error " << lambda - lambda_max
          << std::endl;

  double lambda2;
  u = 1.;

  EigenPower<> mises2(control, mem, -1.5 * lambda);
  mises2.solve(lambda2, A, u);
  deallog << "Eigenvalue " << lambda2 << " Error " << lambda2 - lambda_min
          << std::endl;

  u      = 1.;
  lambda = 0.;
  EigenInverse<> wieland(control, mem);
  wieland.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << " Error " << lambda - lambda_min
          << std::endl;

  u      = 1.;
  lambda = 10.;
  wieland.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << " Error " << lambda - lambda_max
          << std::endl;

  u      = 1.;
  lambda = 10.;
  EigenInverse<> wieland2(control, mem, .2);
  wieland2.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << " Error " << lambda - lambda_max
          << std::endl;
}
