// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that the BatchOptimizer::copy_from() function works as expected.
//
// Here we use only dictionary substitution, and invoke no symbolic
// optimizations.

#include "../tests.h"

#include "sd_common_tests/batch_optimizer_08.h"

int
main()
{
  initlog();
  deallog << std::setprecision(10);

  // Show the difference between a SymEngine "value" and
  // an evaluated, floating point number
  // deallog << std::setprecision(3);

  const enum SD::OptimizerType     opt_method = SD::OptimizerType::dictionary;
  const enum SD::OptimizationFlags opt_flags =
    SD::OptimizationFlags::optimize_all;

  const int dim = 2;
  run_tests<dim, opt_method, opt_flags>();
}
