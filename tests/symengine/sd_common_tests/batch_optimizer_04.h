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

#include <deal.II/differentiation/sd.h>

#include "../../tests.h"

#include "serialization.h"

using namespace dealii;
namespace SD = Differentiation::SD;

template <int dim,
          typename NumberType,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
test_serialization()
{
  namespace SD = Differentiation::SD;
  typedef SD::Expression SD_number_t;

  // Serialization
  // Note: We scope all symbolic variables to make sure that they are
  // totally deleted. They're all persistent pointers in the background,
  // so to fully verify the save-load condition we have to be doubly
  // sure that everything is nuked. We'll lift the basic serialization
  // routing from the serialization header.
  std::ostringstream oss;

  {
    // Define
    const SD_number_t x("x"), y("y");
    const SD_number_t f = x * y;
    const SD_number_t g = x / y;

    // Substitution map
    const SD::types::substitution_map sub_vals =
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

    // Serialize
    deallog << "Serializing..." << std::endl;
    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    oa << optimizer;
  }

  // The output from the serialization of the LLVM optimizer
  // confuses CTest, because it's binary. So don't output it
  // in that cse.
  if (opt_method != SD::OptimizerType::llvm)
    {
      deallog << "Serialized optimizer: " << std::endl;
      deallog << oss.str() << std::endl;
    }
  // const std::string old_optimizer_string = oss.str();

  // Deserialization
  {
    // Define
    const SD_number_t x("x"), y("y");
    const SD_number_t f = x * y;
    const SD_number_t g = x / y;

    // Optimizer
    SD::BatchOptimizer<NumberType> optimizer;

    // Deserialze
    deallog << "Deserializing..." << std::endl;
    std::istringstream            iss(oss.str());
    boost::archive::text_iarchive ia(iss, boost::archive::no_header);
    ia >> optimizer;

    if (opt_method != SD::OptimizerType::llvm)
      {
        // Verify the correctness of the deserialization
        std::ostringstream            oss_2;
        boost::archive::text_oarchive oa(oss_2, boost::archive::no_header);
        oa << optimizer;

        deallog << "Deserialized optimizer: " << std::endl;
        deallog << oss_2.str() << std::endl;
        // const std::string new_optimizer_string = oss_2.str();

        // deallog << "Compare: " <<
        // new_optimizer_string.compare(old_optimizer_string) << std::endl;
        // AssertThrow(!new_optimizer_string.compare(old_optimizer_string),
        // ExcInternalError());
      }

    // Evaluate: Verify the functioning of the deserialized optimizer
    const NumberType val_f = optimizer.evaluate(f);
    const NumberType val_g = optimizer.evaluate(g);
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

  deallog.push("Serialization");
  {
    // deallog.push("Integer");
    // test_tensor<dim,int>(n_runs, timer);
    // deallog.pop();

    deallog.push("Float");
    test_serialization<dim, float, opt_method, opt_flags>();
    deallog.pop();

    deallog.push("Double");
    test_serialization<dim, double, opt_method, opt_flags>();
    deallog.pop();

    // The LLVM optimizer does not currently support complex numbers.
    if (opt_method != SD::OptimizerType::llvm)
      {
        deallog.push("Complex float");
        test_serialization<dim, std::complex<float>, opt_method, opt_flags>();
        deallog.pop();

        deallog.push("Complex double");
        test_serialization<dim, std::complex<double>, opt_method, opt_flags>();
        deallog.pop();
      }
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
