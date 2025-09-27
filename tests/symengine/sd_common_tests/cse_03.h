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


// Header for common subexpression elimination (CSE) test:
// Coupled double-differentiation example using SD-SD

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/differentiation/sd.h>

#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>

#include "utilities.h"

using namespace dealii;
namespace SD = Differentiation::SD;

// Function and its derivatives
template <int dim, typename NumberType>
struct CoupledFunction
{
  // Firstly, here are a whole lot of helper functions to make
  // computing derivatives of the pseudo-energy function easier

  static NumberType
  det_t(const Tensor<2, dim, NumberType> &t)
  {
    return determinant(t);
  }

  static Tensor<2, dim, NumberType>
  ddet_t_dt(const Tensor<2, dim, NumberType> &t)
  {
    return det_t(t) * transpose(invert(t));
  }

  static Tensor<4, dim, NumberType>
  d2det_t_dt_dt(const Tensor<2, dim, NumberType> &t)
  {
    const Tensor<2, dim, NumberType> t_inv = invert(t);
    Tensor<4, dim, NumberType>       dt_inv_trans_dt;
    // https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Tensors#Derivative_of_the_determinant_of_a_tensor
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l)
            dt_inv_trans_dt[i][j][k][l] = -t_inv[l][i] * t_inv[j][k];

    return det_t(t) * outer_product(transpose(t_inv), transpose(t_inv)) +
           det_t(t) * dt_inv_trans_dt;
  }

  static NumberType
  det_t(const SymmetricTensor<2, dim, NumberType> &t)
  {
    return determinant(t);
  }

  static SymmetricTensor<2, dim, NumberType>
  ddet_t_dt(const SymmetricTensor<2, dim, NumberType> &t)
  {
    const SymmetricTensor<2, dim, NumberType> t_inv = invert(t);
    return det_t(t) * t_inv;
  }

  static SymmetricTensor<4, dim, NumberType>
  d2det_t_dt_dt(const SymmetricTensor<2, dim, NumberType> &t)
  {
    const SymmetricTensor<2, dim, NumberType> t_inv = invert(t);
    SymmetricTensor<4, dim, NumberType>       dt_inv_dt;
    // https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Tensors#Derivative_of_the_determinant_of_a_tensor
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = k; l < dim; ++l)
            dt_inv_dt[i][j][k][l] =
              -0.5 * (t_inv[i][k] * t_inv[j][l] + t_inv[i][l] * t_inv[j][k]);

    return det_t(t) * outer_product(t_inv, t_inv) + det_t(t) * dt_inv_dt;
  }

  static NumberType
  v_dot_v(const Tensor<1, dim, NumberType> &v)
  {
    return v * v;
  }

  static Tensor<1, dim, NumberType>
  dv_dot_v_dv(const Tensor<1, dim, NumberType> &v)
  {
    return 2.0 * v;
  }

  static Tensor<2, dim, NumberType>
  d2v_dot_v_dv_dv(const Tensor<1, dim, NumberType> &v)
  {
    static const Tensor<2, dim, NumberType> I(unit_symmetric_tensor<dim>());
    return 2.0 * I;
  }

  // --------------------------------------------------------------------------

  // Here we define a pseudo-energy function that we wish to differentiate twice
  // with respect to all variables

  static const double sf;

  static NumberType
  psi(const Tensor<2, dim, NumberType> &t,
      const Tensor<1, dim, NumberType> &v,
      const NumberType                 &s)
  {
    return std::pow(det_t(t), 2) * std::pow(v_dot_v(v), 3) * std::pow(s, sf);
  };

  static NumberType
  psi_optimised(const NumberType &det_t,
                const NumberType &v_dot_v,
                const NumberType &s)
  {
    return std::pow(det_t, 2) * std::pow(v_dot_v, 3) * std::pow(s, sf);
  };

  // --------------------------------------------------------------------------

  // Here's the hand-calculated expansion of this energy functions first and
  // second derivatives. Lots of fun...

  static Tensor<2, dim, NumberType>
  dpsi_dt(const Tensor<2, dim, NumberType> &t,
          const Tensor<1, dim, NumberType> &v,
          const NumberType                 &s)
  {
    return 2.0 * std::pow(det_t(t), 1) * ddet_t_dt(t) *
           std::pow(v_dot_v(v), 3) * std::pow(s, sf);
  };

  static Tensor<1, dim, NumberType>
  dpsi_dv(const Tensor<2, dim, NumberType> &t,
          const Tensor<1, dim, NumberType> &v,
          const NumberType                 &s)
  {
    return std::pow(det_t(t), 2) * 3.0 * std::pow(v_dot_v(v), 2) *
           dv_dot_v_dv(v) * std::pow(s, sf);
  };

  static NumberType
  dpsi_ds(const Tensor<2, dim, NumberType> &t,
          const Tensor<1, dim, NumberType> &v,
          const NumberType                 &s)
  {
    return std::pow(det_t(t), 2) * std::pow(v_dot_v(v), 3) * sf *
           std::pow(s, sf - 1.0);
  };

  static Tensor<4, dim, NumberType>
  d2psi_dt_dt(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType                 &s)
  {
    return 2.0 * std::pow(v_dot_v(v), 3) *
           (std::pow(det_t(t), 0) * outer_product(ddet_t_dt(t), ddet_t_dt(t)) +
            std::pow(det_t(t), 1) * d2det_t_dt_dt(t)) *
           std::pow(s, sf);
  };

  static Tensor<3, dim, NumberType>
  d2psi_dv_dt(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType                 &s)
  {
    return 2.0 * std::pow(det_t(t), 1) * 3.0 * std::pow(v_dot_v(v), 2) *
           outer_product(ddet_t_dt(t), dv_dot_v_dv(v)) * std::pow(s, sf);
  };

  static Tensor<2, dim, NumberType>
  d2psi_ds_dt(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType                 &s)
  {
    return 2.0 * std::pow(det_t(t), 1) * ddet_t_dt(t) *
           std::pow(v_dot_v(v), 3) * sf * std::pow(s, sf - 1.0);
  };

  static Tensor<3, dim, NumberType>
  d2psi_dt_dv(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType                 &s)
  {
    return 2.0 * std::pow(det_t(t), 1) * 3.0 * std::pow(v_dot_v(v), 2) *
           outer_product(dv_dot_v_dv(v), ddet_t_dt(t)) * std::pow(s, sf);
  };

  static Tensor<2, dim, NumberType>
  d2psi_dv_dv(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType                 &s)
  {
    return std::pow(det_t(t), 2) * 3.0 *
           (2.0 * std::pow(v_dot_v(v), 1) *
              outer_product(dv_dot_v_dv(v), dv_dot_v_dv(v)) +
            std::pow(v_dot_v(v), 2) * d2v_dot_v_dv_dv(v)) *
           std::pow(s, sf);
  };

  static Tensor<1, dim, NumberType>
  d2psi_ds_dv(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType                 &s)
  {
    return std::pow(det_t(t), 2) * 3.0 * std::pow(v_dot_v(v), 2) *
           dv_dot_v_dv(v) * sf * std::pow(s, sf - 1.0);
  };

  static Tensor<2, dim, NumberType>
  d2psi_dt_ds(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType                 &s)
  {
    return 2.0 * std::pow(det_t(t), 1) * ddet_t_dt(t) *
           std::pow(v_dot_v(v), 3) * sf * std::pow(s, sf - 1.0);
  };

  static Tensor<1, dim, NumberType>
  d2psi_dv_ds(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType                 &s)
  {
    return std::pow(det_t(t), 2) * 3.0 * std::pow(v_dot_v(v), 2) *
           dv_dot_v_dv(v) * sf * std::pow(s, sf - 1.0);
  };

  static NumberType
  d2psi_ds_ds(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType                 &s)
  {
    return std::pow(det_t(t), 2) * std::pow(v_dot_v(v), 3) * sf * (sf - 1.0) *
           std::pow(s, sf - 2.0);
  };
};

template <int dim, typename NumberType>
const double CoupledFunction<dim, NumberType>::sf(2.2);

template <int dim,
          typename number_t,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
evaluate_SD_SD_stored_symbols_optimisation(const Tensor<2, dim> &t,
                                           const Tensor<1, dim> &v,
                                           const double         &s)
{
  using SDNumberType = SD::Expression;

  // Compute the symbolic (derivative) tensors once off
  // ahead of time, and store them for each run of the function
  // Symbols
  static const Tensor<2, dim, SDNumberType> symb_t(
    SD::make_tensor_of_symbols<2, dim>("t"));
  static const Tensor<1, dim, SDNumberType> symb_v(
    SD::make_tensor_of_symbols<1, dim>("v"));
  static const SDNumberType symb_s(SD::make_symbol("s"));
  // Value
  static SDNumberType symb_psi;
  // First derivatives
  static Tensor<2, dim, SDNumberType> symb_dpsi_dt;
  static Tensor<1, dim, SDNumberType> symb_dpsi_dv;
  static SDNumberType                 symb_dpsi_ds;
  // Second derivatives
  static Tensor<4, dim, SDNumberType> symb_d2psi_dt_dt;
  static Tensor<3, dim, SDNumberType> symb_d2psi_dv_dt;
  static Tensor<2, dim, SDNumberType> symb_d2psi_ds_dt;
  static Tensor<3, dim, SDNumberType> symb_d2psi_dt_dv;
  static Tensor<2, dim, SDNumberType> symb_d2psi_dv_dv;
  static Tensor<1, dim, SDNumberType> symb_d2psi_ds_dv;
  static Tensor<2, dim, SDNumberType> symb_d2psi_dt_ds;
  static Tensor<1, dim, SDNumberType> symb_d2psi_dv_ds;
  static SDNumberType                 symb_d2psi_ds_ds;

  static SD::BatchOptimizer<double> optimizer(opt_method, opt_flags);

  static bool symbols_initialised = false;
  if (symbols_initialised == false)
    {
      // Function and its derivatives
      using func_sd = CoupledFunction<dim, SDNumberType>;

      symb_psi = func_sd::psi(symb_t, symb_v, symb_s);

      std::cout << std::string(80, '-') << std::endl;
      std::cout << "Dim: " << dim << std::endl;
      print(std::cout, "symb_t", symb_t);
      print(std::cout, "symb_v", symb_v);
      print(std::cout, "symb_s", symb_s);
      print(std::cout, "symb_psi", symb_psi);

      // First derivative
      symb_dpsi_dt = SD::differentiate(symb_psi, symb_t);
      symb_dpsi_dv = SD::differentiate(symb_psi, symb_v);
      symb_dpsi_ds = SD::differentiate(symb_psi, symb_s);
      // Second derivatives
      symb_d2psi_dt_dt = SD::differentiate(symb_dpsi_dt, symb_t);
      symb_d2psi_dv_dt = SD::differentiate(symb_dpsi_dt, symb_v);
      symb_d2psi_ds_dt = SD::differentiate(symb_dpsi_dt, symb_s);
      symb_d2psi_dt_dv = SD::differentiate(symb_dpsi_dv, symb_t);
      symb_d2psi_dv_dv = SD::differentiate(symb_dpsi_dv, symb_v);
      symb_d2psi_ds_dv = SD::differentiate(symb_dpsi_dv, symb_s);
      symb_d2psi_dt_ds = SD::differentiate(symb_dpsi_ds, symb_t);
      symb_d2psi_dv_ds = SD::differentiate(symb_dpsi_ds, symb_v);
      symb_d2psi_ds_ds = SD::differentiate(symb_dpsi_ds, symb_s);

      print(std::cout, "symb_psi", symb_psi);
      print(std::cout, "symb_dpsi_dt", symb_dpsi_dt);
      print(std::cout, "symb_dpsi_dv", symb_dpsi_dv);
      print(std::cout, "symb_dpsi_ds", symb_dpsi_ds);
      print(std::cout, "symb_d2psi_dt_dt", symb_d2psi_dt_dt);
      print(std::cout, "symb_d2psi_dv_dt", symb_d2psi_dv_dt);
      print(std::cout, "symb_d2psi_ds_dt", symb_d2psi_ds_dt);
      print(std::cout, "symb_d2psi_dt_dv", symb_d2psi_dt_dv);
      print(std::cout, "symb_d2psi_dv_dv", symb_d2psi_dv_dv);
      print(std::cout, "symb_d2psi_ds_dv", symb_d2psi_ds_dv);
      print(std::cout, "symb_d2psi_dt_ds", symb_d2psi_dt_ds);
      print(std::cout, "symb_d2psi_dv_ds", symb_d2psi_dv_ds);
      print(std::cout, "symb_d2psi_ds_ds", symb_d2psi_ds_ds);

      const SD::types::substitution_map sub_vals_optim =
        SD::make_symbol_map(symb_t, symb_v, symb_s);
      optimizer.register_symbols(sub_vals_optim);

      // Value
      optimizer.register_function(symb_psi);
      // First derivative
      optimizer.register_function(symb_dpsi_dt);
      optimizer.register_function(symb_dpsi_dv);
      optimizer.register_function(symb_dpsi_ds);
      // Second derivatives
      optimizer.register_function(symb_d2psi_dt_dt);
      optimizer.register_function(symb_d2psi_dv_dt);
      optimizer.register_function(symb_d2psi_ds_dt);
      optimizer.register_function(symb_d2psi_dt_dv);
      optimizer.register_function(symb_d2psi_dv_dv);
      optimizer.register_function(symb_d2psi_ds_dv);
      optimizer.register_function(symb_d2psi_dt_ds);
      optimizer.register_function(symb_d2psi_dv_ds);
      optimizer.register_function(symb_d2psi_ds_ds);

      optimizer.optimize();

      symbols_initialised = true;

      // Unfortunately, we cannot print out the optimizer state. The result
      // is not stable and depends on standard library, optimization
      // parameters and compiler being used. MM '20
      //
      // deallog << "Optimizer" << std::endl;
      // optimizer.print(std::cout.get_file_stream());
    }

  SD::types::substitution_map sub_vals;
  SD::add_to_substitution_map(sub_vals, symb_t, t);
  SD::add_to_substitution_map(sub_vals, symb_v, v);
  SD::add_to_substitution_map(sub_vals, symb_s, s);

  optimizer.substitute(sub_vals);

  // Value
  const double psi = optimizer.evaluate(symb_psi);
  // First derivatives
  const Tensor<2, dim, double> dpsi_dt = optimizer.evaluate(symb_dpsi_dt);
  const Tensor<1, dim, double> dpsi_dv = optimizer.evaluate(symb_dpsi_dv);
  const double                 dpsi_ds = optimizer.evaluate(symb_dpsi_ds);
  // Second derivatives
  const Tensor<4, dim, double> d2psi_dt_dt =
    optimizer.evaluate(symb_d2psi_dt_dt);
  const Tensor<3, dim, double> d2psi_dv_dt =
    optimizer.evaluate(symb_d2psi_dv_dt);
  const Tensor<2, dim, double> d2psi_ds_dt =
    optimizer.evaluate(symb_d2psi_ds_dt);
  const Tensor<3, dim, double> d2psi_dt_dv =
    optimizer.evaluate(symb_d2psi_dt_dv);
  const Tensor<2, dim, double> d2psi_dv_dv =
    optimizer.evaluate(symb_d2psi_dv_dv);
  const Tensor<1, dim, double> d2psi_ds_dv =
    optimizer.evaluate(symb_d2psi_ds_dv);
  const Tensor<2, dim, double> d2psi_dt_ds =
    optimizer.evaluate(symb_d2psi_dt_ds);
  const Tensor<1, dim, double> d2psi_dv_ds =
    optimizer.evaluate(symb_d2psi_dv_ds);
  const double d2psi_ds_ds = optimizer.evaluate(symb_d2psi_ds_ds);

  // Verify the result
  using func              = CoupledFunction<dim, double>;
  static const double tol = 1e-12;
  AssertThrow(std::abs(psi - func::psi(t, v, s)) < tol,
              ExcMessage("No match for function value."));
  AssertThrow(std::abs((dpsi_dt - func::dpsi_dt(t, v, s)).norm()) < tol,
              ExcMessage("No match for first derivative."));
  AssertThrow(std::abs((dpsi_dv - func::dpsi_dv(t, v, s)).norm()) < tol,
              ExcMessage("No match for first derivative."));
  AssertThrow(std::abs(dpsi_ds - func::dpsi_ds(t, v, s)) < tol,
              ExcMessage("No match for first derivative."));
  AssertThrow(std::abs((d2psi_dt_dt - func::d2psi_dt_dt(t, v, s)).norm()) < tol,
              ExcMessage("No match for second derivative."));
  AssertThrow(std::abs((d2psi_dv_dt - func::d2psi_dv_dt(t, v, s)).norm()) < tol,
              ExcMessage("No match for second derivative."));
  AssertThrow(std::abs((d2psi_ds_dt - func::d2psi_ds_dt(t, v, s)).norm()) < tol,
              ExcMessage("No match for second derivative."));
  AssertThrow(std::abs((d2psi_dt_dv - func::d2psi_dt_dv(t, v, s)).norm()) < tol,
              ExcMessage("No match for second derivative."));
  AssertThrow(std::abs((d2psi_dv_dv - func::d2psi_dv_dv(t, v, s)).norm()) < tol,
              ExcMessage("No match for second derivative."));
  AssertThrow(std::abs((d2psi_ds_dv - func::d2psi_ds_dv(t, v, s)).norm()) < tol,
              ExcMessage("No match for second derivative."));
  AssertThrow(std::abs((d2psi_dt_ds - func::d2psi_dt_ds(t, v, s)).norm()) < tol,
              ExcMessage("No match for second derivative."));
  AssertThrow(std::abs((d2psi_dv_ds - func::d2psi_dv_ds(t, v, s)).norm()) < tol,
              ExcMessage("No match for second derivative."));
  AssertThrow(std::abs(d2psi_ds_ds - func::d2psi_ds_ds(t, v, s)) < tol,
              ExcMessage("No match for second derivative."));

  deallog << "OK" << std::endl;
}

template <int                        dim,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
run(const unsigned int n_runs)
{
  deallog.push("Dim " + Utilities::to_string(dim));

  Tensor<2, dim> t = unit_symmetric_tensor<dim>();
  Tensor<1, dim> v;
  double         s = 0.57;
  for (unsigned int i = 0; i < dim; ++i)
    {
      v[i] = 0.275 * (1.0 + i);
      for (unsigned int j = 0; j < dim; ++j)
        t[i][j] = 1.03 + 1.15 * (((i + 1) * j) % dim);
    }
  t *= 0.9;
  v *= 0.63;
  s *= 1.21;

  for (unsigned int i = 0; i < n_runs; ++i)
    evaluate_SD_SD_stored_symbols_optimisation<dim,
                                               double,
                                               opt_method,
                                               opt_flags>(t, v, s);

  deallog.pop();

  deallog << "OK" << std::endl;
}
