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


// Check that serialization for the BatchOptimizer works as expected.
// In this test, we destroy only the optimizer and keep the other expressions
// that were used to initialize it.
//
// Here we invoke the LLVM optimizer before symbolic evaluation takes place,
// and we invoke no additional symbolic optimizations as well.

#include "../tests.h"

#include "sd_common_tests/batch_optimizer_05.h"

int
main()
{
  initlog();
  deallog << std::setprecision(10);

  const enum SD::OptimizerType     opt_method = SD::OptimizerType::llvm;
  const enum SD::OptimizationFlags opt_flags =
    SD::OptimizationFlags::optimize_default;

  const int dim = 2;
  run_tests<dim, opt_method, opt_flags>();
}
