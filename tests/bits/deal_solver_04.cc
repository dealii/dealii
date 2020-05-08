// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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


// test the MINRES solver


#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_qmrs.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <typeinfo>

#include "../tests.h"

#include "../testmatrix.h"

template <typename SolverType,
          typename MatrixType,
          typename VectorType,
          class PRECONDITION>
void
check_solve(SolverType &         solver,
            const SolverControl &solver_control,
            const MatrixType &   A,
            VectorType &         u,
            VectorType &         f,
            const PRECONDITION & P)
{
  deallog << "Solver type: " << typeid(solver).name() << std::endl;

  u = 0.;
  f = 1.;
  try
    {
      solver.solve(A, u, f, P);
    }
  catch (std::exception &e)
    {
      deallog << e.what() << std::endl;
      abort();
    }

  deallog << "Solver stopped after " << solver_control.last_step()
          << " iterations" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  GrowingVectorMemory<> mem;
  SolverControl         control(100, 1.e-3);

  const unsigned int size = 32;
  unsigned int       dim  = (size - 1) * (size - 1);

  deallog << "Size " << size << " Unknowns " << dim << std::endl;

  // Make matrix
  FDMatrix        testproblem(size, size);
  SparsityPattern structure(dim, dim, 5);
  testproblem.five_point_structure(structure);
  structure.compress();
  SparseMatrix<double> A(structure);
  testproblem.five_point(A);

  Vector<double> f(dim);
  Vector<double> u(dim);
  f = 1.;

  SolverMinRes<>       solver(control, mem);
  PreconditionIdentity preconditioner;
  check_solve(solver, control, A, u, f, preconditioner);
}
