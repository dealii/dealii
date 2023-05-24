//-----------------------------------------------------------
//
//    Copyright (C) 2023 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/petsc_ts.templates.h>

#include "../tests.h"

/**
 * Tests user defined Vector and Matrix types
 * and exceptions handling for PETSCWrappers::TimeStepper.
 */

class VectorType : public Subscriptor
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

class MatrixType : public Subscriptor
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

using TimeStepper = PETScWrappers::TimeStepper<VectorType, MatrixType>;

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  VectorType v(nullptr);
  MatrixType A(nullptr);

  TimeStepper myode;

  try
    {
      auto ts = myode.petsc_ts();
      auto t0 = myode.get_time();
      auto dt = myode.get_time_step();
      AssertThrow(ts == static_cast<TS>(myode), ExcInternalError());
      myode.solve(v, A);
    }
  catch (std::exception &exc)
    {
      deallog << exc.what() << std::endl;
    }
}
