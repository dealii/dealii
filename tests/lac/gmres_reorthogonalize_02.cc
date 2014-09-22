// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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


// tests that GMRES builds an orthonormal basis properly for a larger matrix
// and basis size than gmres_orthogonalize_01

#include "../tests.h"
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>



template <typename number>
void test ()
{
  const unsigned int n = 200;
  Vector<number> rhs(n), sol(n);
  rhs = 1.;

  // only add diagonal entries
  SparsityPattern sp(n, n);
  sp.compress();
  SparseMatrix<number> matrix(sp);

  for (unsigned int i=0; i<n; ++i)
    matrix.diag_element(i) = (i+1);

  SolverControl control(1000, 1e2*std::numeric_limits<number>::epsilon());
  typename SolverGMRES<Vector<number> >::AdditionalData data;
  data.max_n_tmp_vectors = 202;

  SolverGMRES<Vector<number> > solver(control, data);
  solver.solve(matrix, sol, rhs, PreconditionIdentity());
}

int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("double");
  test<double>();
  deallog.pop();
  deallog.threshold_double(1.e-4);
  deallog.push("float");
  test<float>();
  deallog.pop();
}

