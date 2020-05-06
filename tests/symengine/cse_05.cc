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


// Another test for the optimizer. Here we compare the output of
// calls to the batch optimizer when using no optimization and when
// invoking CSE.
//
// This test was used to find a subtle bug when the CSE option was used.
// It resulted in a scaling error in the last outputted result / dependent
// variable.

#include <deal.II/differentiation/sd.h>

#include <algorithm>
#include <map>
#include <sstream>
#include <vector>

#include "../tests.h"

using namespace dealii;
namespace SD = dealii::Differentiation::SD;

void
test_derived()
{
  const int            dim = 3;
  const SD::Expression mu_r("mu_r");
  const SD::Expression mu_e("mu_e");

  const SD::Expression kappa("kappa");
  const SD::Expression sf_sat_0("sf_sat_0");
  const SD::Expression H_sat_0("H_sat_0");

  const Tensor<1, dim, SD::Expression> H =
    SD::make_tensor_of_symbols<1, dim>("H");
  const SymmetricTensor<2, dim, SD::Expression> C =
    SD::make_symmetric_tensor_of_symbols<2, dim>("C");

  const SymmetricTensor<2, dim, SD::Expression> C_inv = invert(C);
  const SymmetricTensor<2, dim, SD::Expression> C_bar =
    std::pow(determinant(C), -1.0 / dim) * C;
  const SD::Expression det_C  = determinant(C);
  const SD::Expression J      = std::sqrt(det_C);
  const SD::Expression I1_bar = trace(C_bar);
  const SD::Expression I4     = H * H;
  const SD::Expression I7     = contract3(H, C_inv, H); // H*C_inv*H;

  const SD::Expression psi_vol =
    (kappa / 4.0) * (J * J - 1.0 - 2.0 * std::log(J));
  const SD::Expression psi_sat =
    1.0 +
    sf_sat_0 * std::erf(std::sqrt(numbers::PI) * I4 / (H_sat_0 * H_sat_0));
  const SD::Expression psi_NH = (0.5 * mu_e) * (I1_bar - dim);
  const double         mu_0   = 4.0 * numbers::PI * 1.0e-7;
  const SD::Expression psi_ME = -0.5 * mu_0 * mu_r * J * I7;
  const SD::Expression psi    = psi_vol + psi_sat * psi_NH + psi_ME;

  // These two things should lead to the same result, but when
  // using CSE they do not..
  //  const SD::Expression psi_wrong = H*C*H;
  //  const SD::Expression psi_fixed = contract3(H,C,H);
  //  std::cout << "psi_wrong: " << psi_wrong << std::endl;
  //  std::cout << "psi_fixed: " << psi_fixed << std::endl;
  // // psi_wrong: H_0*(H_0*C_00 + H_1*C_01 + H_2*C_02) + H_1*(H_0*C_01 +
  // H_1*C_11 + H_2*C_12) + H_2*(H_0*C_02 + H_1*C_12 + H_2*C_22)
  // // psi_fixed: H_0**2*C_00 + H_1**2*C_11 + H_2**2*C_22 + 2*H_0*H_1*C_01 +
  // 2*H_0*H_2*C_02 + 2*H_2*H_1*C_12
  //  throw;

  const SymmetricTensor<2, dim, SD::Expression> S = 2.0 * differentiate(psi, C);
  const Tensor<1, dim, SD::Expression>          B = -differentiate(psi, H);

  typename SD::types::substitution_map symbol_value_map;
  SD::add_to_substitution_map(symbol_value_map, mu_r, 5.0);
  SD::add_to_substitution_map(symbol_value_map, mu_e, 30000.0);
  SD::add_to_substitution_map(symbol_value_map, kappa, 14990000.0);
  SD::add_to_substitution_map(symbol_value_map, sf_sat_0, 2.0);
  SD::add_to_substitution_map(symbol_value_map, H_sat_0, 200000.0);
  Tensor<1, dim> H_vals;
  H_vals[2] = 75000;
  //  H_vals[0] = 1000; // New
  //  H_vals[1] = 2000; // New
  SD::add_to_substitution_map(symbol_value_map, H, H_vals);
  SymmetricTensor<2, dim> C_vals;
  C_vals[0][0] = 1.05263;
  C_vals[1][1] = 1.05263;
  C_vals[2][2] = 0.9025;
  //  C_vals[0][1] = 0.05; // New
  //  C_vals[0][2] = 0.1; // New
  //  C_vals[1][2] = 0.15; // New
  SD::add_to_substitution_map(symbol_value_map, C, C_vals);

  auto eval_deal_II =
    [&symbol_value_map, &psi, &S, &B](const bool         with_lambda_opt,
                                      const bool         with_cse,
                                      const unsigned int n_evals) {
      // Configure optimiser
      SD::BatchOptimizer<double> optimiser;
      if (with_lambda_opt == false && with_cse == false)
        optimiser.set_optimization_method(
          SD::OptimizerType::dictionary,
          SD::OptimizationFlags::optimize_default);
      else if (with_lambda_opt == false && with_cse == true)
        optimiser.set_optimization_method(SD::OptimizerType::dictionary,
                                          SD::OptimizationFlags::optimize_cse);
      else if (with_lambda_opt == true && with_cse == false)
        optimiser.set_optimization_method(
          SD::OptimizerType::lambda, SD::OptimizationFlags::optimize_default);
      else // if (with_lambda_opt == true && with_cse == true)
        optimiser.set_optimization_method(SD::OptimizerType::lambda,
                                          SD::OptimizationFlags::optimize_cse);

      // Set independent and dependent variables, then finalise
      optimiser.register_symbols(symbol_value_map);
      //       optimiser.register_function(B);
      optimiser.register_functions(psi, S, B);
      //       optimiser.register_functions(psi,B,S); // Order doesn't influence
      //       result
      optimiser.optimize();

      // Substitute and get result
      std::cout << "Derived: eval_deal_II: "
                << "  with_lambda_opt: " << with_lambda_opt
                << "  with_cse: " << with_cse << "  n_evals: " << n_evals
                << std::endl;
      double res = 0.0;
      for (unsigned int r = 0; r < n_evals; ++r)
        {
          optimiser.substitute(symbol_value_map);
          const std::vector<double> result = optimiser.evaluate();

          for (unsigned int i = 0; i < result.size(); ++i)
            {
              std::cout << "Result[" << i << "]: " << result[i] << std::endl;
              res += result[i];
            }
        }
      deallog << "res: " << std::setprecision(9) << res << std::endl;

      //       std::cout << std::string(100,'-') << std::endl;
      //       optimiser.print(std::cout,with_cse);
      //       std::cout << std::string(100,'-') << std::endl;
    };

  // Perform more than one evaluation to check that all relevant internal
  // data is correctly reset/reinitialized upon the second and later calls
  // to call()/substitute().
  const unsigned int n_evals = 2;

  deallog.push("deal.II: Dict, No CSE");
  {
    eval_deal_II(false, false, n_evals);
  }
  deallog.pop();

  deallog.push("deal.II: Dict, With CSE");
  {
    eval_deal_II(false, true, n_evals);
  }
  deallog.pop();

  deallog.push("deal.II: Lambda, No CSE");
  {
    eval_deal_II(true, false, n_evals);
  }
  deallog.pop();

  deallog.push("deal.II: Lambda, With CSE");
  {
    eval_deal_II(true, true, n_evals);
  }
  deallog.pop();
}

int
main(int argc, char *argv[])
{
  initlog();

  test_derived();
  deallog << "OK" << std::endl;

  return 0;
}
