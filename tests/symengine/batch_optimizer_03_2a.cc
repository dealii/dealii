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


// Check that one can perform substitution of symbolic derivatives that
// are the result of explicit or implicit relationships between symbolic
// variables.
// See tests/symengine/basic_06.cc and tests/symengine/basic_07.cc for
// a more simple example of differentiation of symbols with
// explicit/implicit relationships.
// We invoke the batch optimizer before symbolic evaluation takes place.
//
// Here we use only lambda substitution,and invoke no symbolic optimizations.

#include "../tests.h"

#include "sd_common_tests/batch_optimizer_03.h"

int
main()
{
  initlog();
  deallog << std::setprecision(10);

  const enum SD::OptimizerType     opt_method = SD::OptimizerType::lambda;
  const enum SD::OptimizationFlags opt_flags =
    SD::OptimizationFlags::optimize_default;

  const int dim = 2;
  run_tests<dim, opt_method, opt_flags>();
}
