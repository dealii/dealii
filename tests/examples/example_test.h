// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2023 by the deal.II authors
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

#ifndef dealii_example_test_h
#define dealii_example_test_h

/*
 * This is a stripped down version of the infamous tests.h file suitable
 * for our example step tests. We do not include the much more intrusive
 * full version.
 */

#include <deal.II/base/exceptions.h>
#include <fstream>

/*
 * Given the name of a file, copy it to std::cout
 */

void
cat_file(const char *filename)
{
  {
    std::ifstream in(filename);
    Assert(in, dealii::ExcIO());
    std::cout << in.rdbuf() << "\n";
  }
}


/*
 * Test that a solver converged within a certain range of iteration steps.
 *
 * SolverType_COMMAND is the command to issue, CONTROL_COMMAND a function call
 * that returns the number of iterations (castable to unsigned int), and
 * MIN_ALLOWED, MAX_ALLOWED is the inclusive range of allowed iteration
 * steps.
 */

#define check_solver_within_range(SolverType_COMMAND,                  \
                                  CONTROL_COMMAND,                     \
                                  MIN_ALLOWED,                         \
                                  MAX_ALLOWED)                         \
  {                                                                    \
    try                                                                \
      {                                                                \
        SolverType_COMMAND;                                            \
      }                                                                \
    catch (SolverControl::NoConvergence & exc)                         \
      {}                                                               \
    const unsigned int steps = CONTROL_COMMAND;                        \
    if (steps >= MIN_ALLOWED && steps <= MAX_ALLOWED)                  \
      {                                                                \
        std::cout << "Solver stopped within " << MIN_ALLOWED << " - "  \
                  << MAX_ALLOWED << " iterations" << std::endl;        \
      }                                                                \
    else                                                               \
      {                                                                \
        std::cout << "Solver stopped after " << steps << " iterations" \
                  << std::endl;                                        \
      }                                                                \
  }


#endif // dealii_example_test_h
