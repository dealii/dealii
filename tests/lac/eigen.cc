// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include "testmatrix.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/eigen.h>
#include <deal.II/lac/precondition.h>

int main()
{
  std::ofstream logfile("output");
//  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  GrowingVectorMemory<> mem;
  SolverControl control(1000, 1.e-5);

  const unsigned int size = 10;
  const unsigned int dim = (size-1)*(size-1);

  /*
   * Compute minimal and maximal
   * eigenvalues of the 5-point
   * stencil matrix
   * (Hackbusch:Iterative Lösung...,
   * Satz 4.1.1)
   */
  const double h = 1./size;
  const double s = std::sin(numbers::PI*h/2.);
  const double c = std::cos(numbers::PI*h/2.);
  const double lambda_max = 8.*c*c;
  const double lambda_min = 8.*s*s;

  FDMatrix testproblem(size, size);
  SparsityPattern structure(dim, dim, 5);
  testproblem.five_point_structure(structure);
  structure.compress();
  SparseMatrix<double>  A(structure);
  testproblem.five_point(A);

  Vector<double> u(dim);
  u = 1.;
  double lambda;

  EigenPower<> mises(control, mem);
  mises.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << " Error " << lambda-lambda_max << std::endl;

  double lambda2;
  u = 1.;

  EigenPower<> mises2(control, mem, -1.5*lambda);
  mises2.solve(lambda2, A, u);
  deallog << "Eigenvalue " << lambda2 << " Error " << lambda2-lambda_min << std::endl;

  u = 1.;
  lambda = 0.;
  EigenInverse<> wieland(control, mem);
  wieland.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << " Error " << lambda-lambda_min << std::endl;

  u = 1.;
  lambda = 10.;
  wieland.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << " Error " << lambda-lambda_max << std::endl;

  u = 1.;
  lambda = 10.;
  EigenInverse<> wieland2(control, mem, .2);
  wieland2.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << " Error " << lambda-lambda_max << std::endl;
}
