// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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


// make sure ConsecutiveControl does not allow its (re)usage when not starting
// from 0-th iteration. This test is expected to run, but throw an error.

#include <deal.II/base/utilities.h>

#include <deal.II/lac/solver_control.h>

#include <iostream>

#include "../tests.h"


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4);

  {
    ConsecutiveControl solver_control(12345, 1.e-3, 3, true, true);
    //                             // n_converged_iterations:
    solver_control.check(0, 1.e-1); // 0
    solver_control.check(1, 1.e-4); // 1
    solver_control.check(2, 1.e-5); // 2
    solver_control.check(3, 2.e-3); // 0
    solver_control.check(4, 5.e-4); // 1
    solver_control.check(5, 4.e-4); // 2
    solver_control.check(6, 3.e-4); // 3 ==> success

    // reuse
    solver_control.check(1,
                         3.e-4); // 1  <----- start from 1st iteration NO-NO-NO
    solver_control.check(2, 2.e-4); // 2
    solver_control.check(3, 1.e-4); // 3 ==> success
  }
}
