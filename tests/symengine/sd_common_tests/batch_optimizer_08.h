// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
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

#include <deal.II/differentiation/sd.h>

#include "../../tests.h"

namespace SD = Differentiation::SD;

template <int dim,
          typename NumberType,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
test_copy()
{
  namespace SD      = Differentiation::SD;
  using SD_number_t = SD::Expression;

  // Define
  const SD_number_t x("x"), y("y");
  const SD_number_t f = x * y;
  const SD_number_t g = x / y;

  // Substitution map
  SD::types::substitution_map sub_vals =
    SD::make_substitution_map(std::make_pair(x, NumberType(10.0)),
                              std::make_pair(y, NumberType(2.0)));

  // Optimize
  SD::BatchOptimizer<NumberType> optimizer(opt_method, opt_flags);
  optimizer.register_symbols(sub_vals); // Independent symbols
  optimizer.register_functions(f, g);
  optimizer.optimize();

  // Substitute
  optimizer.substitute(sub_vals);

  // Evaluate
  const NumberType val_f = optimizer.evaluate(f);
  const NumberType val_g = optimizer.evaluate(g);
  deallog << "f: " << f << std::endl
          << "g: " << g << std::endl
          << "f(x,y): " << val_f << std::endl
          << "g(x,y): " << val_g << std::endl
          << std::endl;

  // Define something new
  {
    deallog << "Before copy" << std::endl;

    const SD_number_t x1("x1"), y1("y1");
    const SD_number_t f1 = x1 * y1;
    const SD_number_t g1 = x1 / y1;

    SD::add_to_substitution_map(sub_vals,
                                std::make_pair(x1, NumberType(20.0)),
                                std::make_pair(y1, NumberType(3.0)));

    // Optimizer
    SD::BatchOptimizer<NumberType> new_optimizer;
    new_optimizer.register_symbols(sub_vals); // Independent symbols
    new_optimizer.register_functions(f1, g1);
    new_optimizer.optimize();

    // Substitute
    new_optimizer.substitute(sub_vals);

    // Evaluate
    const NumberType val_f1 = new_optimizer.evaluate(f1);
    const NumberType val_g1 = new_optimizer.evaluate(g1);
    deallog << "f: " << f1 << std::endl
            << "g: " << g1 << std::endl
            << "f(x,y): " << val_f1 << std::endl
            << "g(x,y): " << val_g1 << std::endl
            << std::endl;

    // Perform copy
    deallog << "After copy" << std::endl;
    new_optimizer.copy_from(optimizer);
    Assert(new_optimizer.optimized() == false,
           ExcMessage("Expected new optimizer to be unoptimized."));

    // Extract
    const NumberType val_f = new_optimizer.extract(f, optimizer.evaluate());
    const NumberType val_g = new_optimizer.extract(g, optimizer.evaluate());
    deallog << "f: " << f << std::endl
            << "g: " << g << std::endl
            << "f(x,y): " << val_f << std::endl
            << "g(x,y): " << val_g << std::endl
            << std::endl;
  }

  deallog << "OK" << std::endl;
}


template <int                        dim,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
run_tests()
{
  // Show the difference between a SymEngine "value" and
  // an evaluated, floating point number
  // deallog << std::setprecision(3);

  deallog.push("Copy");
  {
    // deallog.push("Integer");
    // test_tensor<dim,int>(n_runs, timer);
    // deallog.pop();

    deallog.push("Float");
    test_copy<dim, float, opt_method, opt_flags>();
    deallog.pop();

    deallog.push("Double");
    test_copy<dim, double, opt_method, opt_flags>();
    deallog.pop();

    // The LLVM optimizer does not currently support complex numbers.
    if (opt_method != SD::OptimizerType::llvm)
      {
        deallog.push("Complex float");
        test_copy<dim, std::complex<float>, opt_method, opt_flags>();
        deallog.pop();

        deallog.push("Complex double");
        test_copy<dim, std::complex<double>, opt_method, opt_flags>();
        deallog.pop();
      }
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
