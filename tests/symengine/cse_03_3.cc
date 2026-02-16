// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Common subexpression elimination:
// Coupled double-differentiation example using SD-SD (LLVM optimisation)

#include "../tests.h"

#include "sd_common_tests/cse_03.h"

int
main(int argc, char **argv)
{
  initlog();

  const unsigned int n_runs = 2; // Two runs to check for reuse of optimizer
  const enum SD::OptimizerType     opt_method = SD::OptimizerType::llvm;
  const enum SD::OptimizationFlags opt_flags =
    SD::OptimizationFlags::optimize_all;

  run<2, opt_method, opt_flags>(n_runs);
  run<3, opt_method, opt_flags>(n_runs);

  deallog << "OK" << std::endl;
}
