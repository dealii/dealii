// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// Test for symbolic tensor differentiation:
// Coupled double-differentiation example using SD-SD
// This test relates to ad_common_tests/helper_scalar_coupled_4_components_01.h

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

#include "../tests.h"

namespace SD = Differentiation::SD;

// Function and its derivatives
template <int dim, typename NumberType>
struct CoupledFunction
{
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
  v_squ(const Tensor<1, dim, NumberType> &v)
  {
    return v * v;
  }

  static Tensor<1, dim, NumberType>
  dv_squ_dv(const Tensor<1, dim, NumberType> &v)
  {
    return 2.0 * v;
  }

  static Tensor<2, dim, NumberType>
  d2v_squ_dv_dv(const Tensor<1, dim, NumberType> &v)
  {
    static const Tensor<2, dim, NumberType> I(
      unit_symmetric_tensor<dim, NumberType>());
    return 2.0 * I;
  }

  // --------------------------------------------------------------------------

  // Here we define a pseudo-energy function that we wish to differentiate twice
  // with respect to all variables

  static const double sf;

  static NumberType
  psi(const SymmetricTensor<2, dim, NumberType> &st,
      const Tensor<2, dim, NumberType> &         t,
      const Tensor<1, dim, NumberType> &         v,
      const NumberType &                         s)
  {
    return std::pow(det_t(st), 2) * std::pow(det_t(t), 2) *
           std::pow(v_squ(v), 3) * std::pow(s, sf);
  };

  // --------------------------------------------------------------------------

  // Here's the hand-calculated expansion of this energy functions first and
  // second derivatives. Lots of fun...

  static SymmetricTensor<2, dim, NumberType>
  dpsi_dst(const SymmetricTensor<2, dim, NumberType> &st,
           const Tensor<2, dim, NumberType> &         t,
           const Tensor<1, dim, NumberType> &         v,
           const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(st), 1) * ddet_t_dt(st) *
           std::pow(det_t(t), 2) * std::pow(v_squ(v), 3) * std::pow(s, sf);
  };

  static Tensor<2, dim, NumberType>
  dpsi_dt(const SymmetricTensor<2, dim, NumberType> &st,
          const Tensor<2, dim, NumberType> &         t,
          const Tensor<1, dim, NumberType> &         v,
          const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(st), 2) * std::pow(det_t(t), 1) * ddet_t_dt(t) *
           std::pow(v_squ(v), 3) * std::pow(s, sf);
  };

  static Tensor<1, dim, NumberType>
  dpsi_dv(const SymmetricTensor<2, dim, NumberType> &st,
          const Tensor<2, dim, NumberType> &         t,
          const Tensor<1, dim, NumberType> &         v,
          const NumberType &                         s)
  {
    return std::pow(det_t(st), 2) * std::pow(det_t(t), 2) * 3.0 *
           std::pow(v_squ(v), 2) * dv_squ_dv(v) * std::pow(s, sf);
  };

  static NumberType
  dpsi_ds(const SymmetricTensor<2, dim, NumberType> &st,
          const Tensor<2, dim, NumberType> &         t,
          const Tensor<1, dim, NumberType> &         v,
          const NumberType &                         s)
  {
    return std::pow(det_t(st), 2) * std::pow(det_t(t), 2) *
           std::pow(v_squ(v), 3) * sf * std::pow(s, sf - 1.0);
  };

  static SymmetricTensor<4, dim, NumberType>
  d2psi_dst_x_dst(const SymmetricTensor<2, dim, NumberType> &st,
                  const Tensor<2, dim, NumberType> &         t,
                  const Tensor<1, dim, NumberType> &         v,
                  const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(t), 2) * std::pow(v_squ(v), 3) *
           (std::pow(det_t(st), 0) *
              outer_product(ddet_t_dt(st), ddet_t_dt(st)) +
            std::pow(det_t(st), 1) * d2det_t_dt_dt(st)) *
           std::pow(s, sf);
  };

  static Tensor<4, dim, NumberType>
  d2psi_dst_x_dt(const SymmetricTensor<2, dim, NumberType> &st,
                 const Tensor<2, dim, NumberType> &         t,
                 const Tensor<1, dim, NumberType> &         v,
                 const NumberType &                         s)
  {
    return 4.0 * std::pow(det_t(st), 1) * std::pow(det_t(t), 1) *
           outer_product(Tensor<2, dim, NumberType>(ddet_t_dt(st)),
                         ddet_t_dt(t)) *
           std::pow(v_squ(v), 3) * std::pow(s, sf);
  };

  static Tensor<3, dim, NumberType>
  d2psi_dst_x_dv(const SymmetricTensor<2, dim, NumberType> &st,
                 const Tensor<2, dim, NumberType> &         t,
                 const Tensor<1, dim, NumberType> &         v,
                 const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(t), 2) * std::pow(det_t(st), 1) * 3.0 *
           std::pow(v_squ(v), 2) *
           outer_product(Tensor<2, dim, NumberType>(ddet_t_dt(st)),
                         dv_squ_dv(v)) *
           std::pow(s, sf);
  };

  static SymmetricTensor<2, dim, NumberType>
  d2psi_dst_x_ds(const SymmetricTensor<2, dim, NumberType> &st,
                 const Tensor<2, dim, NumberType> &         t,
                 const Tensor<1, dim, NumberType> &         v,
                 const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(st), 1) * ddet_t_dt(st) *
           std::pow(det_t(t), 2) * std::pow(v_squ(v), 3) * sf *
           std::pow(s, sf - 1.0);
  };

  static Tensor<4, dim, NumberType>
  d2psi_dt_x_dst(const SymmetricTensor<2, dim, NumberType> &st,
                 const Tensor<2, dim, NumberType> &         t,
                 const Tensor<1, dim, NumberType> &         v,
                 const NumberType &                         s)
  {
    return 4.0 * std::pow(det_t(st), 1) * std::pow(det_t(t), 1) *
           outer_product(ddet_t_dt(t),
                         Tensor<2, dim, NumberType>(ddet_t_dt(st))) *
           std::pow(v_squ(v), 3) * std::pow(s, sf);
  };

  static Tensor<4, dim, NumberType>
  d2psi_dt_x_dt(const SymmetricTensor<2, dim, NumberType> &st,
                const Tensor<2, dim, NumberType> &         t,
                const Tensor<1, dim, NumberType> &         v,
                const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(st), 2) * std::pow(v_squ(v), 3) *
           (pow(det_t(t), 0) * outer_product(ddet_t_dt(t), ddet_t_dt(t)) +
            std::pow(det_t(t), 1) * d2det_t_dt_dt(t)) *
           std::pow(s, sf);
  };

  static Tensor<3, dim, NumberType>
  d2psi_dt_x_dv(const SymmetricTensor<2, dim, NumberType> &st,
                const Tensor<2, dim, NumberType> &         t,
                const Tensor<1, dim, NumberType> &         v,
                const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(st), 2) * std::pow(det_t(t), 1) * 3.0 *
           std::pow(v_squ(v), 2) * outer_product(ddet_t_dt(t), dv_squ_dv(v)) *
           std::pow(s, sf);
  };

  static Tensor<2, dim, NumberType>
  d2psi_dt_x_ds(const SymmetricTensor<2, dim, NumberType> &st,
                const Tensor<2, dim, NumberType> &         t,
                const Tensor<1, dim, NumberType> &         v,
                const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(st), 2) * std::pow(det_t(t), 1) * ddet_t_dt(t) *
           std::pow(v_squ(v), 3) * sf * std::pow(s, sf - 1.0);
  };

  static Tensor<3, dim, NumberType>
  d2psi_dv_x_dst(const SymmetricTensor<2, dim, NumberType> &st,
                 const Tensor<2, dim, NumberType> &         t,
                 const Tensor<1, dim, NumberType> &         v,
                 const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(t), 2) * std::pow(det_t(st), 1) * 3.0 *
           std::pow(v_squ(v), 2) *
           outer_product(dv_squ_dv(v),
                         Tensor<2, dim, NumberType>(ddet_t_dt(st))) *
           std::pow(s, sf);
  };

  static Tensor<3, dim, NumberType>
  d2psi_dv_x_dt(const SymmetricTensor<2, dim, NumberType> &st,
                const Tensor<2, dim, NumberType> &         t,
                const Tensor<1, dim, NumberType> &         v,
                const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(st), 2) * std::pow(det_t(t), 1) * 3.0 *
           std::pow(v_squ(v), 2) * outer_product(dv_squ_dv(v), ddet_t_dt(t)) *
           std::pow(s, sf);
  };

  static Tensor<2, dim, NumberType>
  d2psi_dv_x_dv(const SymmetricTensor<2, dim, NumberType> &st,
                const Tensor<2, dim, NumberType> &         t,
                const Tensor<1, dim, NumberType> &         v,
                const NumberType &                         s)
  {
    return std::pow(det_t(st), 2) * std::pow(det_t(t), 2) * 3.0 *
           (2.0 * std::pow(v_squ(v), 1) *
              outer_product(dv_squ_dv(v), dv_squ_dv(v)) +
            std::pow(v_squ(v), 2) * d2v_squ_dv_dv(v)) *
           std::pow(s, sf);
  };

  static Tensor<1, dim, NumberType>
  d2psi_dv_x_ds(const SymmetricTensor<2, dim, NumberType> &st,
                const Tensor<2, dim, NumberType> &         t,
                const Tensor<1, dim, NumberType> &         v,
                const NumberType &                         s)
  {
    return std::pow(det_t(st), 2) * std::pow(det_t(t), 2) * 3.0 *
           std::pow(v_squ(v), 2) * dv_squ_dv(v) * sf * std::pow(s, sf - 1.0);
  };

  static SymmetricTensor<2, dim, NumberType>
  d2psi_ds_x_dst(const SymmetricTensor<2, dim, NumberType> &st,
                 const Tensor<2, dim, NumberType> &         t,
                 const Tensor<1, dim, NumberType> &         v,
                 const NumberType &                         s)
  {
    return 2.0 * ddet_t_dt(st) * std::pow(det_t(t), 2) *
           std::pow(det_t(st), 1) * std::pow(v_squ(v), 3) * sf *
           std::pow(s, sf - 1.0);
  };

  static Tensor<2, dim, NumberType>
  d2psi_ds_x_dt(const SymmetricTensor<2, dim, NumberType> &st,
                const Tensor<2, dim, NumberType> &         t,
                const Tensor<1, dim, NumberType> &         v,
                const NumberType &                         s)
  {
    return 2.0 * std::pow(det_t(st), 2) * std::pow(det_t(t), 1) * ddet_t_dt(t) *
           std::pow(v_squ(v), 3) * sf * std::pow(s, sf - 1.0);
  };

  static Tensor<1, dim, NumberType>
  d2psi_ds_x_dv(const SymmetricTensor<2, dim, NumberType> &st,
                const Tensor<2, dim, NumberType> &         t,
                const Tensor<1, dim, NumberType> &         v,
                const NumberType &                         s)
  {
    return std::pow(det_t(st), 2) * std::pow(det_t(t), 2) * 3.0 *
           std::pow(v_squ(v), 2) * dv_squ_dv(v) * sf * std::pow(s, sf - 1.0);
  };

  static NumberType
  d2psi_ds_x_ds(const SymmetricTensor<2, dim, NumberType> &st,
                const Tensor<2, dim, NumberType> &         t,
                const Tensor<1, dim, NumberType> &         v,
                const NumberType &                         s)
  {
    return std::pow(det_t(st), 2) * std::pow(det_t(t), 2) *
           std::pow(v_squ(v), 3) * sf * (sf - 1.0) * std::pow(s, sf - 2.0);
  };
};

template <int dim, typename NumberType>
const double CoupledFunction<dim, NumberType>::sf(2.2);

template <typename Stream, typename NumberType>
void
print(Stream &stream, const std::string &name, const NumberType &val)
{
  stream << name << ": " << val << std::endl;
}
template <typename Stream, int dim, typename NumberType>
void
print(Stream &                          stream,
      const std::string &               name,
      const Tensor<0, dim, NumberType> &val)
{
  stream << name << ": " << val << std::endl;
}
template <typename Stream, int dim, typename NumberType>
void
print(Stream &                          stream,
      const std::string &               name,
      const Tensor<1, dim, NumberType> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    stream << name << "[" << i << "]: " << t[i] << std::endl;
}
template <typename Stream, int dim, typename NumberType>
void
print(Stream &                          stream,
      const std::string &               name,
      const Tensor<2, dim, NumberType> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      stream << name << "[" << i << "][" << j << "]: " << t[i][j] << std::endl;
}
template <typename Stream, int dim, typename NumberType>
void
print(Stream &                          stream,
      const std::string &               name,
      const Tensor<3, dim, NumberType> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        stream << name << "[" << i << "][" << j << "][" << k
               << "]: " << t[i][j][k] << std::endl;
}
template <typename Stream, int dim, typename NumberType>
void
print(Stream &                          stream,
      const std::string &               name,
      const Tensor<4, dim, NumberType> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          stream << name << "[" << i << "][" << j << "][" << k << "][" << l
                 << "]: " << t[i][j][k][l] << std::endl;
}

template <int dim, typename number_t>
void
test_symmetric_tensor_tensor_vector_scalar_coupled(
  const SymmetricTensor<2, dim, number_t> &st,
  const Tensor<2, dim, number_t> &         t,
  const Tensor<1, dim, number_t> &         v,
  const number_t &                         s)
{
  using SDNumberType = SD::Expression;

  const SymmetricTensor<2, dim, SDNumberType> symb_st(
    SD::make_symmetric_tensor_of_symbols<2, dim>("st"));
  const Tensor<2, dim, SDNumberType> symb_t(
    SD::make_tensor_of_symbols<2, dim>("t"));
  const Tensor<1, dim, SDNumberType> symb_v(
    SD::make_vector_of_symbols<dim>("v"));
  const SDNumberType symb_s(SD::make_symbol("s"));

  // Value
  SDNumberType symb_psi;

  // First derivatives
  SymmetricTensor<2, dim, SDNumberType> symb_dpsi_dst;
  Tensor<2, dim, SDNumberType>          symb_dpsi_dt;
  Tensor<1, dim, SDNumberType>          symb_dpsi_dv;
  SDNumberType                          symb_dpsi_ds;

  // Second derivatives
  SymmetricTensor<4, dim, SDNumberType> symb_d2psi_dst_x_dst;
  Tensor<4, dim, SDNumberType>          symb_d2psi_dst_x_dt;
  Tensor<3, dim, SDNumberType>          symb_d2psi_dst_x_dv;
  SymmetricTensor<2, dim, SDNumberType> symb_d2psi_dst_x_ds;

  Tensor<4, dim, SDNumberType> symb_d2psi_dt_x_dst;
  Tensor<4, dim, SDNumberType> symb_d2psi_dt_x_dt;
  Tensor<3, dim, SDNumberType> symb_d2psi_dt_x_dv;
  Tensor<2, dim, SDNumberType> symb_d2psi_dt_x_ds;

  Tensor<3, dim, SDNumberType> symb_d2psi_dv_x_dst;
  Tensor<3, dim, SDNumberType> symb_d2psi_dv_x_dt;
  Tensor<2, dim, SDNumberType> symb_d2psi_dv_x_dv;
  Tensor<1, dim, SDNumberType> symb_d2psi_dv_x_ds;

  SymmetricTensor<2, dim, SDNumberType> symb_d2psi_ds_x_dst;
  Tensor<2, dim, SDNumberType>          symb_d2psi_ds_x_dt;
  Tensor<1, dim, SDNumberType>          symb_d2psi_ds_x_dv;
  SDNumberType                          symb_d2psi_ds_x_ds;

  // Function and its derivatives
  using func_sd = CoupledFunction<dim, SDNumberType>;
  symb_psi      = func_sd::psi(symb_st, symb_t, symb_v, symb_s);

  // First derivative
  symb_dpsi_dst = SD::differentiate(symb_psi, symb_st);
  symb_dpsi_dt  = SD::differentiate(symb_psi, symb_t);
  symb_dpsi_dv  = SD::differentiate(symb_psi, symb_v);
  symb_dpsi_ds  = SD::differentiate(symb_psi, symb_s);

  // Second derivatives
  symb_d2psi_dst_x_dst = SD::differentiate(symb_dpsi_dst, symb_st);
  symb_d2psi_dst_x_dt  = SD::differentiate(symb_dpsi_dst, symb_t);
  symb_d2psi_dst_x_dv  = SD::differentiate(symb_dpsi_dst, symb_v);
  symb_d2psi_dst_x_ds  = SD::differentiate(symb_dpsi_dst, symb_s);

  symb_d2psi_dt_x_dst = SD::differentiate(symb_dpsi_dt, symb_st);
  symb_d2psi_dt_x_dt  = SD::differentiate(symb_dpsi_dt, symb_t);
  symb_d2psi_dt_x_dv  = SD::differentiate(symb_dpsi_dt, symb_v);
  symb_d2psi_dt_x_ds  = SD::differentiate(symb_dpsi_dt, symb_s);

  symb_d2psi_dv_x_dst = SD::differentiate(symb_dpsi_dv, symb_st);
  symb_d2psi_dv_x_dt  = SD::differentiate(symb_dpsi_dv, symb_t);
  symb_d2psi_dv_x_dv  = SD::differentiate(symb_dpsi_dv, symb_v);
  symb_d2psi_dv_x_ds  = SD::differentiate(symb_dpsi_dv, symb_s);

  symb_d2psi_ds_x_dst = SD::differentiate(symb_dpsi_ds, symb_st);
  symb_d2psi_ds_x_dt  = SD::differentiate(symb_dpsi_ds, symb_t);
  symb_d2psi_ds_x_dv  = SD::differentiate(symb_dpsi_ds, symb_v);
  symb_d2psi_ds_x_ds  = SD::differentiate(symb_dpsi_ds, symb_s);

#if DEBUG
  print(std::cout, "symb_st", symb_st);
  print(std::cout, "symb_t", symb_t);
  print(std::cout, "symb_v", symb_v);
  print(std::cout, "symb_s", symb_s);

  print(std::cout, "symb_psi", symb_psi);

  print(std::cout, "symb_dpsi_dst", symb_dpsi_dst);
  print(std::cout, "symb_dpsi_dt", symb_dpsi_dt);
  print(std::cout, "symb_dpsi_dv", symb_dpsi_dv);
  print(std::cout, "symb_dpsi_ds", symb_dpsi_ds);

  print(std::cout, "symb_d2psi_dst_x_dst", symb_d2psi_dst_x_dst);
  print(std::cout, "symb_d2psi_dst_x_dt", symb_d2psi_dst_x_dt);
  print(std::cout, "symb_d2psi_dst_x_dv", symb_d2psi_dst_x_dv);
  print(std::cout, "symb_d2psi_dst_x_ds", symb_d2psi_dst_x_ds);

  print(std::cout, "symb_d2psi_dt_x_dst", symb_d2psi_dt_x_dst);
  print(std::cout, "symb_d2psi_dt_x_dt", symb_d2psi_dt_x_dt);
  print(std::cout, "symb_d2psi_dt_x_dv", symb_d2psi_dt_x_dv);
  print(std::cout, "symb_d2psi_dt_x_ds", symb_d2psi_dt_x_ds);

  print(std::cout, "symb_d2psi_dv_x_dst", symb_d2psi_dv_x_dst);
  print(std::cout, "symb_d2psi_dv_x_dt", symb_d2psi_dv_x_dt);
  print(std::cout, "symb_d2psi_dv_x_dv", symb_d2psi_dv_x_dv);
  print(std::cout, "symb_d2psi_dv_x_ds", symb_d2psi_dv_x_ds);

  print(std::cout, "symb_d2psi_ds_x_dst", symb_d2psi_ds_x_dst);
  print(std::cout, "symb_d2psi_ds_x_dt", symb_d2psi_ds_x_dt);
  print(std::cout, "symb_d2psi_ds_x_dv", symb_d2psi_ds_x_dv);
  print(std::cout, "symb_d2psi_ds_x_ds", symb_d2psi_ds_x_ds);
#endif

  // It pains me to do this manually...
  // This is why we have utility functions to create substitution maps
  // and do evaluation
  typename SD::types::substitution_map sub_vals;
  sub_vals[symb_s] = SDNumberType(s);
  for (unsigned int i = 0; i < dim; ++i)
    {
      sub_vals[symb_v[i]] = SDNumberType(v[i]);
      for (unsigned int j = 0; j < dim; ++j)
        {
          sub_vals[symb_t[i][j]]  = SDNumberType(t[i][j]);
          sub_vals[symb_st[i][j]] = SDNumberType(st[i][j]);
        }
    }

  // Perform substitution and evaluation
  // Value
  number_t psi;

  // First derivatives
  SymmetricTensor<2, dim, number_t> dpsi_dst;
  Tensor<2, dim, number_t>          dpsi_dt;
  Tensor<1, dim, number_t>          dpsi_dv;
  number_t                          dpsi_ds;

  // Second derivatives
  SymmetricTensor<4, dim, number_t> d2psi_dst_x_dst;
  Tensor<4, dim, number_t>          d2psi_dst_x_dt;
  Tensor<3, dim, number_t>          d2psi_dst_x_dv;
  SymmetricTensor<2, dim, number_t> d2psi_dst_x_ds;

  Tensor<4, dim, number_t> d2psi_dt_x_dst;
  Tensor<4, dim, number_t> d2psi_dt_x_dt;
  Tensor<3, dim, number_t> d2psi_dt_x_dv;
  Tensor<2, dim, number_t> d2psi_dt_x_ds;

  Tensor<3, dim, number_t> d2psi_dv_x_dst;
  Tensor<3, dim, number_t> d2psi_dv_x_dt;
  Tensor<2, dim, number_t> d2psi_dv_x_dv;
  Tensor<1, dim, number_t> d2psi_dv_x_ds;

  SymmetricTensor<2, dim, number_t> d2psi_ds_x_dst;
  Tensor<2, dim, number_t>          d2psi_ds_x_dt;
  Tensor<1, dim, number_t>          d2psi_ds_x_dv;
  number_t                          d2psi_ds_x_ds;

  psi     = symb_psi.substitute_and_evaluate<number_t>(sub_vals);
  dpsi_ds = symb_dpsi_ds.substitute_and_evaluate<number_t>(sub_vals);
  d2psi_ds_x_ds =
    symb_d2psi_ds_x_ds.substitute_and_evaluate<number_t>(sub_vals);
  for (unsigned int i = 0; i < dim; ++i)
    {
      dpsi_dv[i] =
        symb_dpsi_dv[i].template substitute_and_evaluate<number_t>(sub_vals);
      d2psi_dv_x_ds[i] =
        symb_d2psi_dv_x_ds[i].template substitute_and_evaluate<number_t>(
          sub_vals);
      d2psi_ds_x_dv[i] =
        symb_d2psi_ds_x_dv[i].template substitute_and_evaluate<number_t>(
          sub_vals);

      for (unsigned int j = 0; j < dim; ++j)
        {
          dpsi_dst[i][j] =
            symb_dpsi_dst[i][j].template substitute_and_evaluate<number_t>(
              sub_vals);
          dpsi_dt[i][j] =
            symb_dpsi_dt[i][j].template substitute_and_evaluate<number_t>(
              sub_vals);
          d2psi_dst_x_ds[i][j] =
            symb_d2psi_dst_x_ds[i][j]
              .template substitute_and_evaluate<number_t>(sub_vals);
          d2psi_dt_x_ds[i][j] =
            symb_d2psi_dt_x_ds[i][j].template substitute_and_evaluate<number_t>(
              sub_vals);
          d2psi_dv_x_dv[i][j] =
            symb_d2psi_dv_x_dv[i][j].template substitute_and_evaluate<number_t>(
              sub_vals);
          d2psi_ds_x_dst[i][j] =
            symb_d2psi_ds_x_dst[i][j]
              .template substitute_and_evaluate<number_t>(sub_vals);
          d2psi_ds_x_dt[i][j] =
            symb_d2psi_ds_x_dt[i][j].template substitute_and_evaluate<number_t>(
              sub_vals);

          for (unsigned int k = 0; k < dim; ++k)
            {
              d2psi_dst_x_dv[i][j][k] =
                symb_d2psi_dst_x_dv[i][j][k]
                  .template substitute_and_evaluate<number_t>(sub_vals);
              d2psi_dt_x_dv[i][j][k] =
                symb_d2psi_dt_x_dv[i][j][k]
                  .template substitute_and_evaluate<number_t>(sub_vals);
              d2psi_dv_x_dst[i][j][k] =
                symb_d2psi_dv_x_dst[i][j][k]
                  .template substitute_and_evaluate<number_t>(sub_vals);
              d2psi_dv_x_dt[i][j][k] =
                symb_d2psi_dv_x_dt[i][j][k]
                  .template substitute_and_evaluate<number_t>(sub_vals);

              for (unsigned int l = 0; l < dim; ++l)
                {
                  d2psi_dst_x_dst[i][j][k][l] =
                    symb_d2psi_dst_x_dst[i][j][k][l]
                      .template substitute_and_evaluate<number_t>(sub_vals);
                  d2psi_dst_x_dt[i][j][k][l] =
                    symb_d2psi_dst_x_dt[i][j][k][l]
                      .template substitute_and_evaluate<number_t>(sub_vals);
                  d2psi_dt_x_dst[i][j][k][l] =
                    symb_d2psi_dt_x_dst[i][j][k][l]
                      .template substitute_and_evaluate<number_t>(sub_vals);
                  d2psi_dt_x_dt[i][j][k][l] =
                    symb_d2psi_dt_x_dt[i][j][k][l]
                      .template substitute_and_evaluate<number_t>(sub_vals);
                }
            }
        }
    }

  // Verify the result
  using func              = CoupledFunction<dim, number_t>;
  static const double tol = 1e-9;
  Assert(std::abs(psi - func::psi(st, t, v, s)) < tol,
         ExcMessage("No match for function value."));

  Assert(std::abs((dpsi_dst - func::dpsi_dst(st, t, v, s)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  Assert(std::abs((dpsi_dt - func::dpsi_dt(st, t, v, s)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  Assert(std::abs((dpsi_dv - func::dpsi_dv(st, t, v, s)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  Assert(std::abs(dpsi_ds - func::dpsi_ds(st, t, v, s)) < tol,
         ExcMessage("No match for first derivative."));

  Assert(std::abs(
           (d2psi_dst_x_dst - func::d2psi_dst_x_dst(st, t, v, s)).norm()) < tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dst_x_dt - func::d2psi_dst_x_dt(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dst_x_dv - func::d2psi_dst_x_dv(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dst_x_ds - func::d2psi_dst_x_ds(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dt_x_dst - func::d2psi_dt_x_dst(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dt_x_dt - func::d2psi_dt_x_dt(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dt_x_dv - func::d2psi_dt_x_dv(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dt_x_ds - func::d2psi_dt_x_ds(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dv_x_dst - func::d2psi_dv_x_dst(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dv_x_dt - func::d2psi_dv_x_dt(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dv_x_dv - func::d2psi_dv_x_dv(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_dv_x_ds - func::d2psi_dv_x_ds(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_ds_x_dst - func::d2psi_ds_x_dst(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_ds_x_dt - func::d2psi_ds_x_dt(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs((d2psi_ds_x_dv - func::d2psi_ds_x_dv(st, t, v, s)).norm()) <
           tol,
         ExcMessage("No match for second derivative."));
  Assert(std::abs(d2psi_ds_x_ds - func::d2psi_ds_x_ds(st, t, v, s)) < tol,
         ExcMessage("No match for second derivative."));
}

template <int dim, typename number_t>
void
run()
{
  deallog.push("dim " + Utilities::to_string(dim));

  SymmetricTensor<2, dim, number_t> st =
    SymmetricTensor<2, dim, number_t>(unit_symmetric_tensor<dim>());
  Tensor<2, dim, number_t> t =
    Tensor<2, dim, number_t>(unit_symmetric_tensor<dim>());
  Tensor<1, dim, number_t> v;
  double                   s = 0.57;
  for (unsigned int i = 0; i < dim; ++i)
    {
      v[i] = 0.275 * (1.0 + i);
      for (unsigned int j = 0; j < dim; ++j)
        {
          t[i][j]  = 1.03 + 1.15 * (((i + 1) * j) % dim);
          st[i][j] = 1.01 + 1.25 * (((i + 2) * j) % dim);
        }
    }
  st *= 0.85;
  t *= 0.9;
  v *= 0.63;
  s *= 1.21;

  test_symmetric_tensor_tensor_vector_scalar_coupled<dim, number_t>(st,
                                                                    t,
                                                                    v,
                                                                    s);

  deallog << "OK" << std::endl;
  deallog.pop();
}

int
main(int argc, char **argv)
{
  initlog();

  deallog.push("Double");
  {
    run<2, double>();
    run<3, double>();
  }
  deallog.pop();

  //  deallog.push("Float");
  //  {
  //    run<2, float>();
  //    run<3, float>();
  //  }
  //  deallog.pop();

  deallog << "OK" << std::endl;
}
