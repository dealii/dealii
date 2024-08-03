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

#ifndef dealii_example_test_h
#define dealii_example_test_h

/*
 * This is a stripped down version of the infamous tests.h file suitable
 * for our example step tests. We do not include the much more intrusive
 * full version.
 */

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iostream>

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

#define check_solver_within_range(                                        \
  OSTREAM, SolverType_COMMAND, CONTROL_COMMAND, MIN_ALLOWED, MAX_ALLOWED) \
  {                                                                       \
    try                                                                   \
      {                                                                   \
        SolverType_COMMAND;                                               \
      }                                                                   \
    catch (SolverControl::NoConvergence &)                                \
      {}                                                                  \
    const unsigned int steps = CONTROL_COMMAND;                           \
    if (steps >= MIN_ALLOWED && steps <= MAX_ALLOWED)                     \
      {                                                                   \
        OSTREAM << "Solver stopped within " << MIN_ALLOWED << " - "       \
                << MAX_ALLOWED << " iterations" << std::endl;             \
      }                                                                   \
    else                                                                  \
      {                                                                   \
        OSTREAM << "Solver stopped after " << steps << " iterations"      \
                << std::endl;                                             \
      }                                                                   \
  }

std::string   deallogname;
std::ofstream deallogfile;

void
initlog(const bool                    console = false,
        const std::ios_base::fmtflags flags   = std::ios::showpoint |
                                              std::ios::left)
{
  deallogname = "output";
  deallogfile.open(deallogname);
  dealii::deallog.attach(deallogfile, true, flags);
  dealii::deallog.depth_console(console ? 10 : 0);
}

#endif // dealii_example_test_h
