// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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


// Check that undefine_macros.h undefines all assertion macros. Check this by
// defining functions with the same names.

#include <deal.II/base/exceptions.h>

#include <deal.II/sundials/arkode.h>

#include "../tests.h"

// must come after any deal.II headers to undefine them
#include <deal.II/base/undefine_macros.h>

int
Assert(int, int)
{
  return 42;
}
int
AssertARKode(int, int)
{
  return 42;
}

int
AssertCuda(int, int)
{
  return 42;
}

int
AssertCusolver(int, int)
{
  return 42;
}

int
AssertCusparse(int, int)
{
  return 42;
}

int
AssertDimension(int, int)
{
  return 42;
}

int
AssertIDA(int, int)
{
  return 42;
}

int
AssertIndexRange(int, int)
{
  return 42;
}

int
AssertIsFinite(int, int)
{
  return 42;
}

int
AssertKINSOL(int, int)
{
  return 42;
}

int
AssertNothrow(int, int)
{
  return 42;
}

int
AssertNothrowCuda(int, int)
{
  return 42;
}

int
AssertNothrowCusparse(int, int)
{
  return 42;
}

int
AssertThrow(int, int)
{
  return 42;
}

int
AssertThrowMPI(int, int)
{
  return 42;
}

int
AssertVectorVectorDimension(int, int)
{
  return 42;
}

int
DeclException0(int, int)
{
  return 42;
}

int
DeclException1(int, int)
{
  return 42;
}

int
DeclException2(int, int)
{
  return 42;
}

int
DeclException3(int, int)
{
  return 42;
}

int
DeclException4(int, int)
{
  return 42;
}

int
DeclException5(int, int)
{
  return 42;
}

int
DeclExceptionMsg(int, int)
{
  return 42;
}


int
main()
{
  initlog();

  deallog << "OK" << std::endl;
}
