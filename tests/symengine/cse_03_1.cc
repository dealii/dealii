// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Common subexpression elimination:
// Coupled double-differentiation example using SD-SD (dictionary substitution)

#include "../tests.h"

#include "sd_common_tests/cse_03.h"

int
main(int argc, char **argv)
{
  initlog();

  const unsigned int n_runs = 2; // Two runs to check for reuse of optimizer
  const enum SD::OptimizerType     opt_method = SD::OptimizerType::dictionary;
  const enum SD::OptimizationFlags opt_flags =
    SD::OptimizationFlags::optimize_all;

  run<2, opt_method, opt_flags>(n_runs);
  run<3, opt_method, opt_flags>(n_runs);

  deallog << "OK" << std::endl;
}
