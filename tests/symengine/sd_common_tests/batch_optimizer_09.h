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


// Check that assertions for some illegal operations when using the
// BatchOptimizer get triggered.

#include <deal.II/differentiation/sd.h>

#include "../../tests.h"

#include "utilities.h"

namespace SD = Differentiation::SD;

template <int dim,
          typename NumberType,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
test_extract_evaluate()
{
  namespace SD           = Differentiation::SD;
  using SD_number_t      = SD::Expression;
  using SD_tensor_t      = Tensor<2, dim, SD_number_t>;
  using SD_symm_tensor_t = SymmetricTensor<2, dim, SD_number_t>;

  // Define
  const SD_number_t      x("x"), y("y");
  const SD_tensor_t      T = SD::make_tensor_of_symbols<2, dim>("T");
  const SD_symm_tensor_t S = SD::make_symmetric_tensor_of_symbols<2, dim>("S");
  const SD_number_t      f = x * y;
  const SD_tensor_t      G = T * (x + y);
  const SD_symm_tensor_t H = S / (x - y);

  // Substitution map
  const SD::types::substitution_map sub_vals = SD::make_substitution_map(
    std::make_pair(x, NumberType(10.0)),
    std::make_pair(y, NumberType(2.0)),
    std::make_pair(T, make_tensor<dim>(NumberType(3.0))),
    std::make_pair(S, make_symm_tensor<dim>(NumberType(4.0))));

  // Optimize
  SD::BatchOptimizer<NumberType> optimizer(opt_method, opt_flags);
  optimizer.register_symbols(sub_vals);  // Independent symbols
  optimizer.register_functions(f, G, H); // Dependent functions
  optimizer.optimize();

  // Substitute
  optimizer.substitute(sub_vals);

  // Evaluate
  deallog.push("Evaluation");
  {
    const NumberType                          val_f = optimizer.evaluate(f);
    const Tensor<2, dim, NumberType>          val_G = optimizer.evaluate(G);
    const SymmetricTensor<2, dim, NumberType> val_H = optimizer.evaluate(H);
    deallog << "f: " << f << std::endl
            << "G: " << G << std::endl
            << "H: " << H << std::endl
            << "f(x,y): " << val_f << std::endl
            << "G(x,y,T): " << val_G << std::endl
            << "H(x,y,S): " << val_H << std::endl
            << std::endl;
  }
  deallog.pop();

  // Now we setup a second optimizer, which will only perform value extraction
  SD::BatchOptimizer<NumberType> extraction_optimizer(opt_method, opt_flags);
  extraction_optimizer.register_functions(f, G, H); // Dependent functions

  // Extract
  deallog.push("Extraction");
  {
    const auto values = optimizer.evaluate();

    const NumberType extr_val_f = extraction_optimizer.extract(f, values);
    const Tensor<2, dim, NumberType> extr_val_G =
      extraction_optimizer.extract(G, values);
    const SymmetricTensor<2, dim, NumberType> extr_val_H =
      extraction_optimizer.extract(H, values);
    deallog << "f(x,y): " << extr_val_f << std::endl
            << "G(x,y,T): " << extr_val_G << std::endl
            << "H(x,y,S): " << extr_val_H << std::endl
            << std::endl;
  }
  deallog.pop();

  // Attempt evaluation with extraction-specific optimiser. We expect this to
  // fail since extraction_optimizer.optimize() hasn't been called.
  deallog.push("Invalid evaluation");
  {
    extraction_optimizer.register_symbols(sub_vals);

    deallog << "Going to try evaluation without substitution..." << std::endl;

    try
      {
        const NumberType val_f = extraction_optimizer.evaluate(f);
        deallog << "f(x,y): " << val_f << std::endl << std::endl;
      }
    catch (std::exception &)
      {
        deallog << "Caught invalid evaluation before substitution: Scalar"
                << std::endl;
      }

    try
      {
        const Tensor<2, dim, NumberType> val_G =
          extraction_optimizer.evaluate(G);
        deallog << "G(x,y,T): " << val_G << std::endl << std::endl;
      }
    catch (std::exception &)
      {
        deallog << "Caught invalid evaluation before substitution: Tensor"
                << std::endl;
      }

    try
      {
        const SymmetricTensor<2, dim, NumberType> val_H =
          extraction_optimizer.evaluate(H);
        deallog << "H(x,y,S): " << val_H << std::endl << std::endl;
      }
    catch (std::exception &)
      {
        deallog
          << "Caught invalid evaluation before substitution: SymmetricTensor"
          << std::endl;
      }

    deallog << "Going to try substitution without optimization..." << std::endl;

    try
      {
        extraction_optimizer.substitute(sub_vals);
      }
    catch (std::exception &)
      {
        deallog << "Caught invalid substitution before optimization."
                << std::endl;
      }

    deallog << "Going to try evaluation without optimization..." << std::endl;

    try
      {
        const NumberType val_f = extraction_optimizer.evaluate(f);
        deallog << "f(x,y): " << val_f << std::endl << std::endl;
      }
    catch (std::exception &)
      {
        deallog << "Caught invalid evaluation before optimization: Scalar"
                << std::endl;
      }

    try
      {
        const Tensor<2, dim, NumberType> val_G =
          extraction_optimizer.evaluate(G);
        deallog << "G(x,y,T): " << val_G << std::endl << std::endl;
      }
    catch (std::exception &)
      {
        deallog << "Caught invalid evaluation before optimization: Tensor"
                << std::endl;
      }

    try
      {
        const SymmetricTensor<2, dim, NumberType> val_H =
          extraction_optimizer.evaluate(H);
        deallog << "H(x,y,S): " << val_H << std::endl << std::endl;
      }
    catch (std::exception &)
      {
        deallog
          << "Caught invalid evaluation before optimization: SymmetricTensor"
          << std::endl;
      }
  }
  deallog.pop();

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

  deallog.push("Double");
  try
    {
      test_extract_evaluate<dim, double, opt_method, opt_flags>();
    }
  catch (...)
    {
      deallog.pop();
      deallog << "NOT OK" << std::endl;
    }
  deallog.pop();

  deallog << "OK" << std::endl;
}
