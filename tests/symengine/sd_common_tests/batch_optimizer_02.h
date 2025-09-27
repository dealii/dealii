// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
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
// batch optimizer before symbolic evaluation takes place.

#include <deal.II/base/timer.h>

#include <deal.II/differentiation/sd.h>

#include "../../tests.h"

#include "utilities.h"

namespace SD = Differentiation::SD;

template <typename NumberType,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
test_scalar(const int n_runs, TimerOutput &timer)
{
  std::cout << std::string(80, '-') << std::endl;
  std::cout << "Scalar" << std::endl;

  using SD_number_t = SD::Expression;

  const NumberType  a = NumberType(1.5);
  const SD_number_t x("x");
  const SD_number_t y("y");
  timer.enter_subsection("Value calculation");
  const SD_number_t symb_s = a + NumberType(2.0) * std::pow(x, y) +
                             a * y / std::log(x) + std::sin(x * y * y / a);
  std::cout << "symb_s: " << symb_s << std::endl;
  timer.leave_subsection("Value calculation");

  timer.enter_subsection("Differentiation");
  const SD_number_t symb_ds_dx = SD::differentiate(symb_s, x);
  std::cout << "symb_ds_dx: " << symb_ds_dx << std::endl;
  timer.leave_subsection("Differentiation");

  SD::types::substitution_map sub_vals;
  SD::add_to_substitution_map(sub_vals, SD::make_substitution_map(x, 2.5));
  SD::add_to_substitution_map(sub_vals, SD::make_substitution_map(y, 1.5));

  deallog.push("Substitution");
  {
    TimerOutput::Scope timer_scope(timer, "Unoptimised substitution");
    for (unsigned int i = 0; i < n_runs; ++i)
      {
        const SD_number_t subs_s     = SD::substitute(symb_s, sub_vals);
        const SD_number_t subs_ds_dx = SD::substitute(symb_ds_dx, sub_vals);
        if (i == 0)
          std::cout << "substitute: "
                    << "  s: " << subs_s << "  ds_dx: " << subs_ds_dx
                    << std::endl;

        const NumberType val_s     = static_cast<NumberType>(subs_s);
        const NumberType val_ds_dx = static_cast<NumberType>(subs_ds_dx);
        if (i == 0)
          std::cout << "evaluation: "
                    << "  s: " << val_s << "  ds_dx: " << val_ds_dx
                    << std::endl;
      }
  }
  deallog.pop();

  deallog.push("Optimisation + substitution");
  {
    // Send our symbolic expression through
    // a batch optimiser
    timer.enter_subsection("Optimisation");
    SD::BatchOptimizer<NumberType> optimizer(opt_method, opt_flags);
    optimizer.register_symbols(sub_vals); // Independent symbols
    optimizer.register_function(
      symb_s); // Dependent symbolic expressions (scalar)
    optimizer.register_function(
      symb_ds_dx); // Dependent symbolic expressions (scalar)
    optimizer.optimize();
    timer.leave_subsection("Optimisation");

    TimerOutput::Scope timer_scope(timer, "Optimised substitution");
    for (unsigned int i = 0; i < n_runs; ++i)
      {
        optimizer.substitute(sub_vals);

        const NumberType val_s     = optimizer.evaluate(symb_s);
        const NumberType val_ds_dx = optimizer.evaluate(symb_ds_dx);
        if (i == 0)
          std::cout << "evaluation: "
                    << "  s: " << val_s << "  ds_dx: " << val_ds_dx
                    << std::endl;

        // The result is not stable and depends on standard library, number
        // type, optimization parameters and compiler being used. However, we
        // also don't want to maintain a whole bunch of blesses output
        // variants. So we set a loose tolerance here.
        static const double tol = 1e-3;
        {
          const NumberType blessed_val_s     = 11.2896853345;
          const NumberType blessed_val_ds_dx = 2.44062401526;

          AssertThrow(std::abs(val_s - blessed_val_s) < tol,
                      ExcMessage("No match for function value."));
          AssertThrow(std::abs(val_ds_dx - blessed_val_ds_dx) < tol,
                      ExcMessage("No match for first derivative."));
        }
      }
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}

template <int dim,
          typename NumberType,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
test_tensor(const int n_runs, TimerOutput &timer)
{
  std::cout << std::string(80, '-') << std::endl;
  std::cout << "Tensor: Dim: " << dim << std::endl;

  namespace SD           = Differentiation::SD;
  using SD_number_t      = SD::Expression;
  using SD_tensor_t      = Tensor<2, dim, SD_number_t>;
  using SD_symm_tensor_t = SymmetricTensor<2, dim, SD_number_t>;

  const NumberType       a = NumberType(1.5);
  const SD_number_t      x("x");
  const SD_tensor_t      y(SD::make_tensor_of_symbols<2, dim>("y"));
  const SD_symm_tensor_t z(SD::make_symmetric_tensor_of_symbols<2, dim>("z"));

  timer.enter_subsection("Value calculation");
  const SD_number_t symb_s =
    a + NumberType(2.0) * std::pow(x, determinant(y)) +
    a * determinant(y) * std::log((z * z) / determinant(y)) +
    std::sin((z * symmetrize(y)) / a);
  std::cout << "symb_s: " << symb_s << std::endl;
  timer.leave_subsection("Value calculation");

  timer.enter_subsection("Differentiation");
  const SD_number_t      symb_ds_dx = SD::differentiate(symb_s, x);
  const SD_tensor_t      symb_ds_dy = SD::differentiate(symb_s, y);
  const SD_symm_tensor_t symb_ds_dz = SD::differentiate(symb_s, z);
  std::cout << "symb_ds_dy: " << symb_ds_dy << std::endl;
  timer.leave_subsection("Differentiation");

  SD::types::substitution_map sub_vals;
  SD::add_to_substitution_map(sub_vals, SD::make_substitution_map(x, 2.5));
  SD::add_to_substitution_map(
    sub_vals, SD::make_substitution_map(y, make_tensor<dim>(NumberType(2.2))));
  SD::add_to_substitution_map(
    sub_vals,
    SD::make_substitution_map(z, make_symm_tensor<dim>(NumberType(3.7))));

  deallog.push("Substitution");
  {
    TimerOutput::Scope timer_scope(timer, "Unoptimised substitution");
    for (unsigned int i = 0; i < n_runs; ++i)
      {
        const SD_number_t subs_s     = SD::substitute(symb_s, sub_vals);
        const SD_number_t subs_ds_dx = SD::substitute(symb_ds_dx, sub_vals);
        const SD_tensor_t subs_ds_dy = SD::substitute(symb_ds_dy, sub_vals);
        const SD_symm_tensor_t subs_ds_dz =
          SD::substitute(symb_ds_dz, sub_vals);
        if (i == 0)
          std::cout << "substitute: "
                    << "  s: " << subs_s << "  ds_dx: " << subs_ds_dx
                    << "  ds_dy: " << subs_ds_dy << "  ds_dz: " << subs_ds_dz
                    << std::endl;

        const NumberType val_s     = static_cast<NumberType>(subs_s);
        const NumberType val_ds_dx = static_cast<NumberType>(subs_ds_dx);
        const Tensor<2, dim, NumberType> val_ds_dy =
          static_cast<Tensor<2, dim, NumberType>>(subs_ds_dy);
        const SymmetricTensor<2, dim, NumberType> val_ds_dz =
          static_cast<SymmetricTensor<2, dim, NumberType>>(subs_ds_dz);
        if (i == 0)
          std::cout << "evaluation: "
                    << "  s: " << val_s << "  ds_dx: " << val_ds_dx
                    << "  ds_dy: " << val_ds_dy << "  ds_dz: " << val_ds_dz
                    << std::endl;

        // The result is not stable and depends on standard library, number
        // type, optimization parameters and compiler being used. However, we
        // also don't want to maintain a whole bunch of blesses output
        // variants. So we set a loose tolerance here.
        static const double tol = 1e-3;
        if (dim == 2)
          {
            const NumberType                 blessed_val_s     = 733.145160811;
            const NumberType                 blessed_val_ds_dx = 1803.37488238;
            const Tensor<2, dim, NumberType> blessed_val_ds_dy(
              {{4021.22068318, -2721.50003932},
               {-2073.13387957, 2076.40954282}});
            const SymmetricTensor<2, dim, NumberType> blessed_val_ds_dz =
              symmetrize(
                Tensor<2, dim, NumberType>({{1.67425190103, 1.97851016739},
                                            {1.97851016739, 2.69685737073}}));

            AssertThrow(std::abs(val_s - blessed_val_s) < tol,
                        ExcMessage("No match for function value."));
            AssertThrow(std::abs(val_ds_dx - blessed_val_ds_dx) < tol,
                        ExcMessage("No match for first derivative."));
            AssertThrow((val_ds_dy - blessed_val_ds_dy).norm() < tol,
                        ExcMessage("No match for first derivative."));
            AssertThrow((val_ds_dz - blessed_val_ds_dz).norm() < tol,
                        ExcMessage("No match for first derivative."));
          }
        else
          {
            AssertThrow(false, ExcNotImplemented());
          }
      }
  }
  deallog.pop();

  deallog.push("Optimisation + substitution");
  {
    // Send our symbolic expression through
    // a batch optimiser
    timer.enter_subsection("Optimisation");
    SD::BatchOptimizer<NumberType> optimizer(opt_method, opt_flags);
    optimizer.register_symbols(sub_vals); // Independent symbols
    optimizer.register_function(
      symb_s); // Dependent symbolic expressions (scalar)
    optimizer.register_function(
      symb_ds_dx); // Dependent symbolic expressions (scalar)
    optimizer.register_function(
      symb_ds_dy); // Dependent symbolic expressions (tensor)
    optimizer.register_function(
      symb_ds_dz); // Dependent symbolic expressions (symm tensor)
    optimizer.optimize();

    timer.leave_subsection("Optimisation");


    TimerOutput::Scope timer_scope(timer, "Optimised substitution");
    for (unsigned int i = 0; i < n_runs; ++i)
      {
        optimizer.substitute(sub_vals);

        const NumberType val_s     = optimizer.evaluate(symb_s);
        const NumberType val_ds_dx = optimizer.evaluate(symb_ds_dx);
        const Tensor<2, dim, NumberType> val_ds_dy =
          optimizer.evaluate(symb_ds_dy);
        const SymmetricTensor<2, dim, NumberType> val_ds_dz =
          optimizer.evaluate(symb_ds_dz);

        if (i == 0)
          std::cout << "evaluation: "
                    << "  s: " << val_s << "  ds_dx: " << val_ds_dx
                    << "  ds_dy: " << val_ds_dy << "  ds_dz: " << val_ds_dz
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
run_tests(const int n_runs)
{
  // Show the difference between a SymEngine "value" and
  // an evaluated, floating point number
  // deallog << std::setprecision(12);

  deallog.push("Scalar");
  {
    std::cout << "Scalar values" << std::endl;

    TimerOutput timer(std::cout, TimerOutput::summary, TimerOutput::wall_times);

    // deallog.push("Integer");
    // test_scalar<int>(n_runs, timer);
    // deallog.pop();

    deallog.push("Float");
    test_scalar<float, opt_method, opt_flags>(n_runs, timer);
    deallog.pop();
    //
    deallog.push("Double");
    test_scalar<double, opt_method, opt_flags>(n_runs, timer);
    deallog.pop();

    // Not available yet
    // SymEngine::SymEngineException: Invalid Format: Expected Integer or
    // Rational deallog.push("Complex double");
    // test_tensor<dim,std::complex<double>>();
  }
  deallog.pop();

  deallog.push("Tensor");
  {
    std::cout << "Tensor values" << std::endl;

    TimerOutput timer(std::cout, TimerOutput::summary, TimerOutput::wall_times);

    // deallog.push("Integer");
    // test_tensor<dim,int>(n_runs, timer);
    // deallog.pop();

    deallog.push("Float");
    test_tensor<dim, float, opt_method, opt_flags>(n_runs, timer);
    deallog.pop();

    deallog.push("Double");
    test_tensor<dim, double, opt_method, opt_flags>(n_runs, timer);
    deallog.pop();

    // Not available yet
    // SymEngine::SymEngineException: Invalid Format: Expected Integer or
    // Rational deallog.push("Complex double");
    // test_tensor<dim,std::complex<double>>();
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
