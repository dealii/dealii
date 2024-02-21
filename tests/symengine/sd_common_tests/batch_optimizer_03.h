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


// Check that one can perform substitution of symbolic derivatives that
// are the result of explicit or implicit relationships between symbolic
// variables.
// See tests/symengine/basic_06.cc and tests/symengine/basic_07.cc for
// a more simple example of differentiation of symbols with
// explicit/implicit relationships.
// We invoke the batch optimizer before symbolic evaluation takes place.

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/sd.h>

#include <iostream>
#include <string>

#include "../../tests.h"

#include "utilities.h"

namespace SD = Differentiation::SD;

namespace Diff_Test
{
  template <int dim,
            typename NumberType,
            enum SD::OptimizerType     opt_method,
            enum SD::OptimizationFlags opt_flags>
  void
  test_standard()
  {
    SD::BatchOptimizer<NumberType> optimizer;

    const SD::Expression mu_e_SD(SD::make_symbol("mu_e"));
    const SD::Expression kappa_e_SD(SD::make_symbol("kappa_e"));
    const SymmetricTensor<2, dim, SD::Expression> C_SD(
      SD::make_symmetric_tensor_of_symbols<2, dim>("C"));

    const SD::Expression symbolic_psi =
      0.5 * mu_e_SD * (trace(C_SD) - 3) +
      kappa_e_SD * (std::log(determinant(C_SD)));

    // Compute partial derivatives of energy function
    const SymmetricTensor<2, dim, SD::Expression> symbolic_S =
      2.0 * SD::differentiate(symbolic_psi, C_SD); // S = 2*dpsi_dC
    const SymmetricTensor<4, dim, SD::Expression> symbolic_HH =
      2.0 * SD::differentiate(symbolic_S, C_SD); // HH = 2*dS_dC = 4*d2psi_dC_dC

    // Configure optimizer
    SD::types::substitution_map sub_vals_optim;
    SD::add_to_symbol_map(sub_vals_optim, mu_e_SD, kappa_e_SD, C_SD);

    optimizer.set_optimization_method(opt_method, opt_flags);

    optimizer.register_symbols(sub_vals_optim); // Independent symbols
    optimizer.register_functions(symbolic_psi,
                                 symbolic_S,
                                 symbolic_HH); // Dependent symbolic functions
    optimizer.optimize();

    const bool print_symbols = false;
    if (print_symbols == true)
      {
        print(std::cout, "psi", symbolic_psi);
        print(std::cout, "S", symbolic_S);
        print(std::cout, "HH", symbolic_HH);
      }

    // Numerical substitution
    SD::types::substitution_map               sub_vals_unresolved;
    const double                              mu_e    = 3;
    const double                              kappa_e = 10;
    const SymmetricTensor<2, dim, NumberType> C =
      2.0 * unit_symmetric_tensor<dim>();
    SD::add_to_substitution_map(sub_vals_unresolved,
                                std::make_pair(mu_e_SD, mu_e),
                                std::make_pair(kappa_e_SD, kappa_e),
                                std::make_pair(C_SD, C));
    // NOTE: The recursive substitution is not really required in this case, but
    // good to use in practise in case a more complex energy function is
    // employed later
    const SD::types::substitution_map sub_vals =
      SD::resolve_explicit_dependencies(sub_vals_unresolved);

    // Perform substitution of symbols
    optimizer.substitute(sub_vals);

    // Extract values
    const bool print_values = true;
    if (print_values == true)
      {
        print(deallog, "psi", optimizer.evaluate(symbolic_psi));
        print(deallog, "S", optimizer.evaluate(symbolic_S));
        print(deallog, "HH", optimizer.evaluate(symbolic_HH));
      }
  }

  template <int dim,
            typename NumberType,
            enum SD::OptimizerType     opt_method,
            enum SD::OptimizationFlags opt_flags>
  void
  test_explicit()
  {
    SD::BatchOptimizer<NumberType> optimizer;

    // Set up the primary independent variable
    const SD::Expression mu_e_SD(SD::make_symbol("mu_e"));
    const SD::Expression kappa_e_SD(SD::make_symbol("kappa_e"));
    const SymmetricTensor<2, dim, SD::Expression> C_SD(
      SD::make_symmetric_tensor_of_symbols<2, dim>("C"));

    // Set up an internal variable that treated as an independent variable
    const SD::Expression mu_v_SD(SD::make_symbol("mu_v"));
    const SD::Expression kappa_v_SD(SD::make_symbol("kappa_v"));
    const SymmetricTensor<2, dim, SD::Expression> Qi_SD(
      SD::make_symmetric_tensor_of_symbols<2, dim>("Qi"));

    // Compute the dependent variable
    const SD::Expression symbolic_psi_CQi =
      0.5 * mu_e_SD * (trace(C_SD) - 3) +
      kappa_e_SD * (std::log(determinant(C_SD))) +
      0.5 * mu_v_SD * ((C_SD * Qi_SD) - 3) +
      kappa_v_SD * (std::log(C_SD * Qi_SD));

    // Compute partial derivatives of energy function
    const SymmetricTensor<2, dim, SD::Expression> symbolic_S_CQi =
      2.0 * SD::differentiate(symbolic_psi_CQi, C_SD); // S = 2*dpsi_dC

    // Here is the explicit definition of Q in terms of C
    // We would know this when Q is linearly dependent on C
    SD::types::substitution_map sub_vals_explicit;
    SD::add_to_substitution_map(sub_vals_explicit,
                                std::make_pair(Qi_SD, 2 * C_SD));

    // Note: After performing these operations we should not longer
    // need to consider the internal variable Q at all
    const SD::Expression symbolic_psi =
      SD::substitute(symbolic_psi_CQi, sub_vals_explicit);
    const SymmetricTensor<2, dim, SD::Expression> symbolic_S =
      SD::substitute(symbolic_S_CQi, sub_vals_explicit);
    const SymmetricTensor<4, dim, SD::Expression> symbolic_HH =
      2.0 * SD::differentiate(symbolic_S, C_SD); // HH = 2*dS_dC = 4*d2psi_dC_dC

    // Configure optimizer
    SD::types::substitution_map sub_vals_optim;
    SD::add_to_symbol_map(
      sub_vals_optim, mu_e_SD, kappa_e_SD, C_SD, mu_v_SD, kappa_v_SD);

    optimizer.set_optimization_method(opt_method, opt_flags);

    optimizer.register_symbols(sub_vals_optim); // Independent symbols
    optimizer.register_functions(symbolic_psi,
                                 symbolic_S,
                                 symbolic_HH); // Dependent symbolic functions
    optimizer.optimize();

    const bool print_symbols = false;
    if (print_symbols == true)
      {
        print(std::cout, "psi (Q)", symbolic_psi_CQi);
        print(std::cout, "psi", symbolic_psi);
        print(std::cout, "S (Q)", symbolic_S_CQi);
        print(std::cout, "S", symbolic_S);
        print(std::cout, "HH", symbolic_HH);
      }

    // Numerical substitution
    SD::types::substitution_map               sub_vals_unresolved;
    const double                              mu_e    = 3;
    const double                              kappa_e = 10;
    const SymmetricTensor<2, dim, NumberType> C =
      2.0 * unit_symmetric_tensor<dim>();
    const double mu_v    = 5;
    const double kappa_v = 2;
    SD::add_to_substitution_map(sub_vals_unresolved,
                                std::make_pair(mu_e_SD, mu_e),
                                std::make_pair(kappa_e_SD, kappa_e),
                                std::make_pair(C_SD, C),
                                std::make_pair(mu_v_SD, mu_v),
                                std::make_pair(kappa_v_SD, kappa_v));
    // NOTE: The recursive substitution is not really required in this case, but
    // good to use in practise in case a more complex energy function is
    // employed later
    const SD::types::substitution_map sub_vals =
      SD::resolve_explicit_dependencies(sub_vals_unresolved);

    // Perform substitution of symbols
    optimizer.substitute(sub_vals);

    // Extract values
    const bool print_values = true;
    if (print_values == true)
      {
        print(deallog, "psi", optimizer.evaluate(symbolic_psi));
        print(deallog, "S", optimizer.evaluate(symbolic_S));
        print(deallog, "HH", optimizer.evaluate(symbolic_HH));
      }
  }

  template <int dim,
            typename NumberType,
            enum SD::OptimizerType     opt_method,
            enum SD::OptimizationFlags opt_flags>
  void
  test_implicit()
  {
    SD::BatchOptimizer<NumberType> optimizer;

    // Set up the primary independent variable
    const SD::Expression mu_e_SD(SD::make_symbol("mu_e"));
    const SD::Expression kappa_e_SD(SD::make_symbol("kappa_e"));
    const SymmetricTensor<2, dim, SD::Expression> C_SD(
      SD::make_symmetric_tensor_of_symbols<2, dim>("C"));

    // Set up an internal variable that treated as an independent variable
    const SD::Expression mu_v_SD(SD::make_symbol("mu_v"));
    const SD::Expression kappa_v_SD(SD::make_symbol("kappa_v"));
    const SymmetricTensor<2, dim, SD::Expression> Qi_SD(
      SD::make_symmetric_tensor_of_symbols<2, dim>("Qi"));

    // Compute the dependent variable
    const SD::Expression symbolic_psi =
      0.5 * mu_e_SD * (trace(C_SD) - 3) +
      kappa_e_SD * (std::log(determinant(C_SD))) +
      0.5 * mu_v_SD * ((C_SD * Qi_SD) - 3) +
      kappa_v_SD * (std::log(C_SD * Qi_SD));

    // Compute partial derivatives of energy function
    const SymmetricTensor<2, dim, SD::Expression> symbolic_S =
      2.0 * SD::differentiate(symbolic_psi, C_SD); // S = 2*dpsi_dC

    // Here is the implicit definition of Q in terms of C
    // We would have to define the relationship in this manner when
    // Q is nonlinearly dependent on C
    const SD::types::substitution_map sub_vals_func_symb_dependencies =
      SD::make_symbol_map(C_SD);
    const SymmetricTensor<2, dim, SD::Expression> Qd_SD =
      SD::make_symmetric_tensor_of_symbolic_functions<2, dim>(
        "Qd", sub_vals_func_symb_dependencies);
    const SymmetricTensor<4, dim, SD::Expression> dQd_SD_dC =
      SD::differentiate(Qd_SD, C_SD);
    //    print(deallog,"Qd_SD",Qd_SD);
    //    print(deallog,"dQd_SD_dC",dQd_SD_dC);

    // Now we substitute out the independent internal variable
    // for one that has a sensitivity on the primary independent variable
    const SD::types::substitution_map sub_vals_implicit =
      SD::make_substitution_map(std::make_pair(Qi_SD, Qd_SD));
    const SymmetricTensor<2, dim, SD::Expression> symbolic_S_subs_Q =
      SD::substitute(symbolic_S, sub_vals_implicit);
    // And differentiate. This computes the total derivative of S
    // with respect to C, but accommodates the implicit relationship
    // between Q(C) and C that is yet to be resolved.
    const SymmetricTensor<4, dim, SD::Expression> symbolic_HH_total_impl =
      2.0 * SD::differentiate(symbolic_S_subs_Q,
                              C_SD); // HH = 2*dS_dC = 4*d2psi_dC_dC
    //    print(deallog,"symbolic_HH_total_impl",symbolic_HH_total_impl);

    // An example of an implicit relationship would be
    // dQ/dt = 2.0*C*log( det(C)*log(Q(C) )
    // Since Q(C) is an implicit relationship, we cannot
    // symbolically resolve Q in terms of the components of C:
    // this would be done numerically using a nonlinear
    // solution scheme. That same nonlinear solver
    // can numerically compute dQ_dC for us, so we will assume
    // that we have both of these sets of values.
    const SymmetricTensor<4, dim, SD::Expression> dQ_dC_i_SD =
      SD::make_symmetric_tensor_of_symbols<4, dim>("dQ_dC_i");
    const SD::types::substitution_map sub_vals_explicit =
      SD::make_substitution_map(std::make_pair(Qd_SD, Qi_SD),
                                std::make_pair(SD::differentiate(Qd_SD, C_SD),
                                               dQ_dC_i_SD));
    //    print(deallog,"dQ_dC_i",dQ_dC_i);

    // Substitute derivative terms
    SD::resolve_explicit_dependencies(sub_vals_explicit);
    const SymmetricTensor<4, dim, SD::Expression> symbolic_HH =
      SD::substitute(symbolic_HH_total_impl, sub_vals_explicit);

    // Configure optimizer
    const SD::types::substitution_map sub_vals_optim = SD::make_symbol_map(
      mu_e_SD, kappa_e_SD, C_SD, mu_v_SD, kappa_v_SD, Qi_SD, dQ_dC_i_SD);

    optimizer.set_optimization_method(opt_method, opt_flags);

    optimizer.register_symbols(sub_vals_optim); // Independent symbols
    optimizer.register_functions(symbolic_psi,
                                 symbolic_S,
                                 symbolic_HH); // Dependent symbolic functions
    optimizer.optimize();

    const bool print_symbols = false;
    if (print_symbols == true)
      {
        print(std::cout, "psi", symbolic_psi);
        print(std::cout, "S", symbolic_S);
        print(std::cout, "HH", symbolic_HH);
      }

    // Numerical substitution
    const double                              mu_e    = 3;
    const double                              kappa_e = 10;
    const SymmetricTensor<2, dim, NumberType> C =
      2.0 * unit_symmetric_tensor<dim>();
    const double                              mu_v    = 5;
    const double                              kappa_v = 2;
    const SymmetricTensor<2, dim, NumberType> Q =
      1.5 * unit_symmetric_tensor<dim>();
    const SymmetricTensor<4, dim, NumberType> dQ_dC =
      3.0 * identity_tensor<dim>();
    const SD::types::substitution_map sub_vals_unresolved =
      SD::make_substitution_map(std::make_pair(mu_e_SD, mu_e),
                                std::make_pair(kappa_e_SD, kappa_e),
                                std::make_pair(C_SD, C),
                                std::make_pair(mu_v_SD, mu_v),
                                std::make_pair(kappa_v_SD, kappa_v),
                                std::make_pair(Qi_SD, Q),
                                std::make_pair(dQ_dC_i_SD, dQ_dC));
    // NOTE: The recursive substitution is not really required in this case, but
    // good to use in practise in case a more complex energy function is
    // employed later
    const SD::types::substitution_map sub_vals =
      SD::resolve_explicit_dependencies(sub_vals_unresolved);

    // Perform substitution of symbols
    optimizer.substitute(sub_vals);

    // Extract values
    const bool print_values = true;
    if (print_values == true)
      {
        print(deallog, "psi", optimizer.evaluate(symbolic_psi));
        print(deallog, "S", optimizer.evaluate(symbolic_S));
        print(deallog, "HH", optimizer.evaluate(symbolic_HH));
      }
  }

} // namespace Diff_Test

template <int                        dim,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
run_tests()
{
  using namespace Diff_Test;

  deallog.push("Standard");
  test_standard<dim, double, opt_method, opt_flags>();
  deallog.pop();

  std::cout << std::string(50, '-') << std::endl;
  deallog << std::string(50, '-') << std::endl;

  deallog.push("Explicit");
  test_explicit<dim, double, opt_method, opt_flags>();
  deallog.pop();

  std::cout << std::string(50, '-') << std::endl;
  deallog << std::string(50, '-') << std::endl;

  deallog.push("Implicit");
  test_implicit<dim, double, opt_method, opt_flags>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
