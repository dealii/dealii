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


// Test that all of the low-level math operations function as expected
// and that their values and derivatives can be computed using the
// optimizer.
// For this test the optimizer uses plain old dictionary substitution,
// but all other optimizations (e.g. CSE) have been enabled.
// This test is based on function_verification_01.cc

#include "../tests.h"

#include "sd_common_tests/batch_optimizer_01.h"


int
main()
{
  initlog();
  deallog << std::setprecision(7);

  const enum SD::OptimizerType     opt_method = SD::OptimizerType::dictionary;
  const enum SD::OptimizationFlags opt_flags =
    SD::OptimizationFlags::optimize_all;

  run_tests<opt_method, opt_flags>();

  deallog << "OK" << std::endl;
}
