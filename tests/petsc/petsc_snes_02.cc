// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/petsc_snes.templates.h>

#include "../tests.h"

/**
 * Tests user defined Vector and Matrix types
 * and exceptions handling for PETSCWrappers::NonlinearSolver.
 */

class VectorType : public EnableObserverPointer
{
public:
  explicit VectorType(Vec v)
    : v(v)
  {}

  Vec &
  petsc_vector()
  {
    return v;
  }

private:
  Vec v;
};

class MatrixType : public EnableObserverPointer
{
public:
  explicit MatrixType(Mat A)
    : A(A)
  {}

  Mat &
  petsc_matrix()
  {
    return A;
  }

private:
  Mat A;
};

using NonlinearSolver = PETScWrappers::NonlinearSolver<VectorType, MatrixType>;

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  VectorType v(nullptr);
  MatrixType A(nullptr);

  NonlinearSolver mysolver;

  try
    {
      auto snes = mysolver.petsc_snes();
      AssertThrow(snes == static_cast<SNES>(mysolver), ExcInternalError());
      mysolver.solve(v, A);
    }
  catch (const StandardExceptions::ExcFunctionNotProvided &)
    {
      deallog << "catching expected exception" << std::endl;
    }
  catch (const std::exception &exc)
    {
      deallog << exc.what() << std::endl;
    }
}
