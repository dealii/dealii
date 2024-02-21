// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that the wrapper for symengine numbers can be integrated into the
// tensor class and works as expected.
// This test is the same as symengine_wrapper_03.cc, except that we invoke the
// batch optimizer  before symbolic evaluation takes place.
//
// Here we use only dictionary substitution, but invoke all symbolic
// optimizations.

#include "../tests.h"

#include "sd_common_tests/batch_optimizer_02.h"

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

  const int dim    = 2;
  const int n_runs = 1000;
  run_tests<dim, opt_method, opt_flags>(n_runs);
}
