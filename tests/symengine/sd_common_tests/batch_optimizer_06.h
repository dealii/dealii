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


// Check that serialization for the BatchOptimizer works as expected.
// In this test, we destroy the optimizer as well as all other expressions
// and symbols that were used in the creation of the optimizer.
//
// This test is based off of symengine/sd_common_tests/batch_optimizer_05.h

#include <deal.II/base/numbers.h>
#include <deal.II/base/timer.h>

#include <deal.II/differentiation/sd.h>

#include "../../tests.h"

#include "serialization.h"
#include "utilities.h"


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

  // Serialization
  // Note: We scope all symbolic variables to make sure that they are
  // totally deleted. They're all persistent pointers in the background,
  // so to fully verify the save-load condition we have to be doubly
  // sure that everything is nuked. We'll lift the basic serialization
  // routing from the serialization header.
  std::ostringstream oss;

#define CREATE_SYMBOLS_FUNCTIONS_SUBSTITUTION_MAP                              \
  const NumberType       a = NumberType(1.5);                                  \
  const SD_number_t      x("x");                                               \
  const SD_tensor_t      y(SD::make_tensor_of_symbols<2, dim>("y"));           \
  const SD_symm_tensor_t z(SD::make_symmetric_tensor_of_symbols<2, dim>("z")); \
                                                                               \
  timer.enter_subsection("Value calculation");                                 \
  const SD_number_t symb_s =                                                   \
    a + NumberType(2.0) * std::pow(x, determinant(y)) +                        \
    a * determinant(y) * std::log((z * z) / determinant(y)) +                  \
    std::sin((z * symmetrize(y)) / a);                                         \
  std::cout << "symb_s: " << symb_s << std::endl;                              \
  timer.leave_subsection("Value calculation");                                 \
                                                                               \
  timer.enter_subsection("Differentiation");                                   \
  const SD_number_t      symb_ds_dx = SD::differentiate(symb_s, x);            \
  const SD_tensor_t      symb_ds_dy = SD::differentiate(symb_s, y);            \
  const SD_symm_tensor_t symb_ds_dz = SD::differentiate(symb_s, z);            \
  std::cout << "symb_ds_dy: " << symb_ds_dy << std::endl;                      \
  timer.leave_subsection("Differentiation");                                   \
                                                                               \
  SD::types::substitution_map sub_vals;                                        \
  SD::add_to_substitution_map(sub_vals, SD::make_substitution_map(x, 2.5));    \
  SD::add_to_substitution_map(                                                 \
    sub_vals,                                                                  \
    SD::make_substitution_map(y, make_tensor<dim>(NumberType(2.2))));          \
  SD::add_to_substitution_map(                                                 \
    sub_vals,                                                                  \
    SD::make_substitution_map(z, make_symm_tensor<dim>(NumberType(3.7))));

  // We'll store some results and compare the outcome
  // of the same calculations after deserialization.
  NumberType                          val_s     = NumberType(0);
  NumberType                          val_ds_dx = NumberType(0);
  Tensor<2, dim, NumberType>          val_ds_dy;
  SymmetricTensor<2, dim, NumberType> val_ds_dz;

  // Serialize
  {
    // Do all of the preliminaries
    CREATE_SYMBOLS_FUNCTIONS_SUBSTITUTION_MAP

    deallog.push("Optimisation + substitution");
    {
      // Send our symbolic expression through
      // a batch optimiser
      timer.enter_subsection("Optimisation");
      SD::BatchOptimizer<NumberType> optimizer(opt_method, opt_flags);
      optimizer.register_symbols(sub_vals); // Independent symbols
      optimizer.register_functions(
        symb_s,
        symb_ds_dx,
        symb_ds_dy,
        symb_ds_dz); // Dependent symbolic expressions
      optimizer.optimize();
      timer.leave_subsection("Optimisation");

      timer.enter_subsection("Optimised substitution");
      for (unsigned int i = 0; i < n_runs; ++i)
        {
          optimizer.substitute(sub_vals);

          // Evaluate
          val_s     = optimizer.evaluate(symb_s);
          val_ds_dx = optimizer.evaluate(symb_ds_dx);
          val_ds_dy = optimizer.evaluate(symb_ds_dy);
          val_ds_dz = optimizer.evaluate(symb_ds_dz);

          if (i == 0)
            {
              std::cout << "evaluation: "
                        << "  s: " << val_s << "  ds_dx: " << val_ds_dx
                        << "  ds_dy: " << val_ds_dy << "  ds_dz: " << val_ds_dz
                        << std::endl;
            }
        }
      timer.leave_subsection("Optimised substitution");

      // Evaluate
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

      // Serialize
      std::cout << "Serializing..." << std::endl;
      boost::archive::text_oarchive oa(oss, boost::archive::no_header);
      oa << optimizer;
    }

    // The output from the serialization of the LLVM optimizer
    // confuses CTest, because it's binary. So don't output it
    // in that cse.
    if (opt_method != SD::OptimizerType::llvm)
      {
        std::cout << "Serialized optimizer: " << std::endl;
        std::cout << oss.str() << std::endl;
      }

    // timer.enter_subsection("Serialization");
    {
      // Do all of the preliminaries
      CREATE_SYMBOLS_FUNCTIONS_SUBSTITUTION_MAP

      // Optimizer
      SD::BatchOptimizer<NumberType> optimizer;

      // Deserialize
      std::cout << "Deserializing..." << std::endl;
      std::istringstream            iss(oss.str());
      boost::archive::text_iarchive ia(iss, boost::archive::no_header);
      ia >> optimizer;

      if (opt_method != SD::OptimizerType::llvm)
        {
          // Verify the correctness of the deserialization
          std::ostringstream            oss_2;
          boost::archive::text_oarchive oa(oss_2, boost::archive::no_header);
          oa << optimizer;

          std::cout << "Deserialized optimizer: " << std::endl;
          std::cout << oss_2.str() << std::endl;
        }

      std::cout << "Checking deserialisation..." << std::endl;
      // Check that the original settings persist.
      Assert(optimizer.optimization_flags() == opt_flags, ExcInternalError());
      if (!(opt_method == SD::OptimizerType::llvm &&
            numbers::NumberTraits<NumberType>::is_complex))
        {
          // We cannot do this check under all circumstances, as internally the
          // optimizer type switches.
          Assert(optimizer.optimization_method() == opt_method,
                 ExcInternalError());
        }

      {
        // Check that new optimizer still produces correct results
        // directly from evaluation
        const NumberType new_val_s     = optimizer.evaluate(symb_s);
        const NumberType new_val_ds_dx = optimizer.evaluate(symb_ds_dx);
        const Tensor<2, dim, NumberType> new_val_ds_dy =
          optimizer.evaluate(symb_ds_dy);
        const SymmetricTensor<2, dim, NumberType> new_val_ds_dz =
          optimizer.evaluate(symb_ds_dz);

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
        optimizer.substitute(sub_vals);

        const NumberType new_val_s     = optimizer.evaluate(symb_s);
        const NumberType new_val_ds_dx = optimizer.evaluate(symb_ds_dx);
        const Tensor<2, dim, NumberType> new_val_ds_dy =
          optimizer.evaluate(symb_ds_dy);
        const SymmetricTensor<2, dim, NumberType> new_val_ds_dz =
          optimizer.evaluate(symb_ds_dz);

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
    // timer.leave_subsection("Serialization");
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
