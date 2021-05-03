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


// Check that serialization for the BatchOptimizer works as expected.
// In this test, we destroy only the optimizer and keep the other expressions
// that were used to initialize it.
//
// This test is based off of symengine/sd_common_tests/batch_optimizer_02.h

#include <deal.II/base/timer.h>

#include <deal.II/differentiation/sd.h>

#include "../../tests.h"

#include "serialization.h"
#include "utilities.h"

using namespace dealii;
namespace SD = Differentiation::SD;

template <int dim,
          typename NumberType,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
test_serialization(const int n_runs, TimerOutput &timer)
{
  std::cout << std::string(80, '-') << std::endl;
  std::cout << "Dim: " << dim << std::endl;

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

  deallog.push("Optimisation + substitution");
  {
    // Send our symbolic expression through
    // a batch optimiser
    timer.enter_subsection("Optimisation");
    SD::BatchOptimizer<NumberType> optimizer(opt_method, opt_flags);
    optimizer.register_symbols(sub_vals); // Independent symbols
    optimizer.register_functions(symb_s,
                                 symb_ds_dx,
                                 symb_ds_dy,
                                 symb_ds_dz); // Dependent symbolic expressions
    optimizer.optimize();
    timer.leave_subsection("Optimisation");

    timer.enter_subsection("Optimised substitution");
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
          {
            std::cout << "evaluation: "
                      << "  s: " << val_s << "  ds_dx: " << val_ds_dx
                      << "  ds_dy: " << val_ds_dy << "  ds_dz: " << val_ds_dz
                      << std::endl;
          }
      }
    timer.leave_subsection("Optimised substitution");

    timer.enter_subsection("Serialization");
    {
      const NumberType val_s     = optimizer.evaluate(symb_s);
      const NumberType val_ds_dx = optimizer.evaluate(symb_ds_dx);
      const Tensor<2, dim, NumberType> val_ds_dy =
        optimizer.evaluate(symb_ds_dy);
      const SymmetricTensor<2, dim, NumberType> val_ds_dz =
        optimizer.evaluate(symb_ds_dz);

      std::cout << "Evaluation (pre-serialization): "
                << "  s: " << val_s << "  ds_dx: " << val_ds_dx
                << "  ds_dy: " << val_ds_dy << "  ds_dz: " << val_ds_dz
                << std::endl;

      // The result is not stable and depends on standard library, number
      // type, optimization parameters and compiler being used. However, we
      // also don't want to maintain a whole bunch of blesses output
      // variants. So we set a loose tolerance here.
      constexpr double tol = 1e-2;
      if (dim == 2)
        {
          const NumberType                 blessed_val_s     = 733.145160811;
          const NumberType                 blessed_val_ds_dx = 1803.37488238;
          const Tensor<2, dim, NumberType> blessed_val_ds_dy(
            {{4021.22068318, -2721.50003932}, {-2073.13387957, 2076.40954282}});
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

      // From serialization.h
      std::cout << "Serializing..." << std::endl;
      SD::BatchOptimizer<NumberType> new_optimizer;
      if (opt_method == SD::OptimizerType::llvm)
        {
          verify_no_logging(optimizer, new_optimizer);
        }
      else
        {
          // verify(optimizer, new_optimizer); // Output not stable
          verify_no_logging(optimizer, new_optimizer);
        }

      std::cout << "Checking deserialisation..." << std::endl;
      {
        // Check that the original settings persist.
        Assert(new_optimizer.optimization_method() ==
                 optimizer.optimization_method(),
               ExcInternalError());
        Assert(new_optimizer.optimization_flags() ==
                 optimizer.optimization_flags(),
               ExcInternalError());

        // Check that new optimizer still produces correct results
        // directly from evaluation
        const NumberType new_val_s     = new_optimizer.evaluate(symb_s);
        const NumberType new_val_ds_dx = new_optimizer.evaluate(symb_ds_dx);
        const Tensor<2, dim, NumberType> new_val_ds_dy =
          new_optimizer.evaluate(symb_ds_dy);
        const SymmetricTensor<2, dim, NumberType> new_val_ds_dz =
          new_optimizer.evaluate(symb_ds_dz);

        std::cout << "Evaluation (post-serialization): "
                  << "  s: " << new_val_s << "  ds_dx: " << new_val_ds_dx
                  << "  ds_dy: " << new_val_ds_dy
                  << "  ds_dz: " << new_val_ds_dz << std::endl;

        constexpr double tol = 1e-9;
        AssertThrow(std::abs(new_val_s - val_s) < tol,
                    ExcMessage(
                      "Problem with optimizer function: Serialization"));
        AssertThrow(std::abs(new_val_ds_dx - val_ds_dx) < tol,
                    ExcMessage(
                      "Problem with optimizer function: Serialization"));
        AssertThrow((new_val_ds_dy - val_ds_dy).norm() < tol,
                    ExcMessage(
                      "Problem with optimizer function: Serialization"));
        AssertThrow((new_val_ds_dz - val_ds_dz).norm() < tol,
                    ExcMessage(
                      "Problem with optimizer function: Serialization"));
      }

      // Check that new optimizer still produces correct results
      // directly after substitution (we'll take a shortcut and
      // use the old substitution map).
      {
        new_optimizer.substitute(sub_vals);

        const NumberType new_val_s     = new_optimizer.evaluate(symb_s);
        const NumberType new_val_ds_dx = new_optimizer.evaluate(symb_ds_dx);
        const Tensor<2, dim, NumberType> new_val_ds_dy =
          new_optimizer.evaluate(symb_ds_dy);
        const SymmetricTensor<2, dim, NumberType> new_val_ds_dz =
          new_optimizer.evaluate(symb_ds_dz);

        constexpr double tol = 1e-9;
        AssertThrow(std::abs(new_val_s - val_s) < tol,
                    ExcMessage(
                      "Problem with optimizer function: Serialization"));
        AssertThrow(std::abs(new_val_ds_dx - val_ds_dx) < tol,
                    ExcMessage(
                      "Problem with optimizer function: Serialization"));
        AssertThrow((new_val_ds_dy - val_ds_dy).norm() < tol,
                    ExcMessage(
                      "Problem with optimizer function: Serialization"));
        AssertThrow((new_val_ds_dz - val_ds_dz).norm() < tol,
                    ExcMessage(
                      "Problem with optimizer function: Serialization"));
      }
    }
    timer.leave_subsection("Serialization");
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}


template <int                        dim,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
run_tests(const int n_runs = 1)
{
  // Show the difference between a SymEngine "value" and
  // an evaluated, floating point number
  // deallog << std::setprecision(12);

  deallog.push("Serialization");
  {
    TimerOutput timer(std::cout, TimerOutput::summary, TimerOutput::wall_times);

    // deallog.push("Integer");
    // test_tensor<dim,int>(n_runs, timer);
    // deallog.pop();

    deallog.push("Float");
    test_serialization<dim, float, opt_method, opt_flags>(n_runs, timer);
    deallog.pop();

    deallog.push("Double");
    test_serialization<dim, double, opt_method, opt_flags>(n_runs, timer);
    deallog.pop();

    // The LLVM optimizer does not currently support complex numbers.
    if (opt_method != SD::OptimizerType::llvm)
      {
        deallog.push("Complex float");
        test_serialization<dim, std::complex<float>, opt_method, opt_flags>(
          n_runs, timer);
        deallog.pop();

        deallog.push("Complex double");
        test_serialization<dim, std::complex<double>, opt_method, opt_flags>(
          n_runs, timer);
        deallog.pop();
      }
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
