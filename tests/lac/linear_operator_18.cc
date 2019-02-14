// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

// Test that we cannot accidentally create a linear operator from a
// temporary object

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"

#define CATCH(call)                                           \
  try                                                         \
    {                                                         \
      call;                                                   \
      deallog << "Error: Something went wrong!" << std::endl; \
    }                                                         \
  catch (ExceptionBase & e)                                   \
    {                                                         \
      deallog << e.get_exc_name() << std::endl;               \
    }

using namespace dealii;

int
main()
{
  initlog();
  deal_II_exceptions::disable_abort_on_exception();

  SparseMatrix<double> A;
  const auto           op_A = linear_operator(A);

  CATCH(linear_operator(SparseMatrix<double>()))
  CATCH(linear_operator(SparseMatrix<double>(), A))
  CATCH(linear_operator(A, SparseMatrix<double>()))
  CATCH(linear_operator(SparseMatrix<double>(), SparseMatrix<double>()))

  SolverControl            solver_control;
  SolverCG<Vector<double>> solver(solver_control);
  CATCH(inverse_operator(op_A, solver, SparseMatrix<double>()))
}
