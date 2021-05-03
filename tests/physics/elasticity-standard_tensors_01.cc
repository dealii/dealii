// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// test standard tensor definitions

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/physics/elasticity/standard_tensors.h>

#include "../tests.h"

using namespace dealii::Physics::Elasticity;

const double c10 = 10.0;
const double c01 = 20.0;

// Consider a Mooney-Rivlin material:
// psi = c10.(I1 - dim) + c01.(I2 - dim)
// where I1 = I1(C), I2 = I2(C)
// Then dI1/dC = I and dI2/dC = I1.I - C
// and S = 2 dpsi/dC
template <int dim>
SymmetricTensor<2, dim>
get_S(const Tensor<2, dim> &F)
{
  const SymmetricTensor<2, dim> C  = symmetrize(transpose(F) * F);
  const double                  I1 = first_invariant(C);
  return 2.0 * c10 * StandardTensors<dim>::I +
         2.0 * c01 * (I1 * StandardTensors<dim>::I - C);
}
// For isotropic media, tau = 2.b.dpsi/db == 2.dpsi/db . b
// where b = F.F^{T}, I1(b) == I1(C) and
//  I2(b) == I2(C)
template <int dim>
SymmetricTensor<2, dim>
get_tau(const Tensor<2, dim> &F)
{
  const SymmetricTensor<2, dim> b  = symmetrize(F * transpose(F));
  const double                  I1 = first_invariant(b);
  const SymmetricTensor<2, dim> tmp =
    2.0 * c10 * StandardTensors<dim>::I +
    2.0 * c01 * (I1 * StandardTensors<dim>::I - b);
  return symmetrize(static_cast<Tensor<2, dim>>(tmp) *
                    static_cast<Tensor<2, dim>>(b));
}

template <int dim>
void
test_standard_tensors()
{
  SymmetricTensor<2, dim> t;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      t[i][j] = dim * i + j + 1.0;

  // Check second-order identity tensor I:
  AssertThrow(std::fabs(StandardTensors<dim>::I * t - trace(t)) < 1e-14,
              ExcInternalError());
  AssertThrow(std::fabs(t * StandardTensors<dim>::I - trace(t)) < 1e-14,
              ExcInternalError());

  // Check fourth-order identity tensor II:
  AssertThrow(std::fabs((StandardTensors<dim>::S * t - t).norm()) < 1e-14,
              ExcInternalError());
  AssertThrow(std::fabs((t * StandardTensors<dim>::S - t).norm()) < 1e-14,
              ExcInternalError());

  // Check fourth-order tensor IxI:
  AssertThrow(std::fabs((StandardTensors<dim>::IxI * t -
                         trace(t) * unit_symmetric_tensor<dim>())
                          .norm()) < 1e-14,
              ExcInternalError());
  AssertThrow(std::fabs((t * StandardTensors<dim>::IxI -
                         trace(t) * unit_symmetric_tensor<dim>())
                          .norm()) < 1e-14,
              ExcInternalError());

  // Check spatial deviatoric tensor dev_P:
  AssertThrow(std::fabs((StandardTensors<dim>::dev_P * t -
                         (t - (trace(t) / dim) * unit_symmetric_tensor<dim>()))
                          .norm()) < 1e-14,
              ExcInternalError());
  AssertThrow(std::fabs((t * StandardTensors<dim>::dev_P -
                         (t - (trace(t) / dim) * unit_symmetric_tensor<dim>()))
                          .norm()) < 1e-14,
              ExcInternalError());

  // Check referential deviatoric tensor Dev_P:
  Tensor<2, dim> F(unit_symmetric_tensor<dim>());
  F[0][1]                    = 0.5;
  F[1][0]                    = 0.25;
  const Tensor<2, dim> F_inv = invert(F);

  // Pull-back a fictitious stress tensor, project it onto a deviatoric space,
  // and then push it forward again tau = F.S.F^{T} --> S = F^{-1}*tau*F^{-T}
  const SymmetricTensor<2, dim> s =
    symmetrize(F_inv * static_cast<Tensor<2, dim>>(t) * transpose(F_inv));
  const SymmetricTensor<2, dim> Dev_P_T_x_s =
    StandardTensors<dim>::Dev_P_T(F) * s;
  const SymmetricTensor<2, dim> s_x_Dev_P = s * StandardTensors<dim>::Dev_P(F);

  // Note: The extra factor J^{2/dim} arises due to the definition of Dev_P
  //       including the factor J^{-2/dim}. Ultimately the stress definitions
  //       for s,t do not align with those required to have the direct
  //       relationship s*Dev_P == dev_P*t. For this we would need S
  //       = 2.dW/dC|_{C=\bar{C}} and t = F.S.F^{T} and \bar{C} =
  //       det(F)^{-2/dim} F^{T}.F .
  AssertThrow(std::fabs((symmetrize(std::pow(determinant(F), 2.0 / dim) * F *
                                    static_cast<Tensor<2, dim>>(s_x_Dev_P) *
                                    transpose(F)) -
                         StandardTensors<dim>::dev_P * t)
                          .norm()) < 1e-14,
              ExcInternalError());
  AssertThrow(std::fabs((symmetrize(std::pow(determinant(F), 2.0 / dim) * F *
                                    static_cast<Tensor<2, dim>>(Dev_P_T_x_s) *
                                    transpose(F)) -
                         StandardTensors<dim>::dev_P * t)
                          .norm()) < 1e-14,
              ExcInternalError());

  // Repeat the above exercise for a "real" material response
  const Tensor<2, dim> F_bar = std::pow(determinant(F), -1.0 / dim) * F;
  const SymmetricTensor<2, dim> S_bar = get_S(F_bar);
  const SymmetricTensor<2, dim> tau_bar =
    symmetrize(F_bar * static_cast<Tensor<2, dim>>(S_bar) *
               transpose(F_bar)); // Note: tau_bar = tau(F) |_{F = F_bar}
  AssertThrow(std::fabs((tau_bar - get_tau(F_bar)).norm()) < 1e-9,
              ExcInternalError());
  const SymmetricTensor<2, dim> S_iso = S_bar * StandardTensors<dim>::Dev_P(F);
  const SymmetricTensor<2, dim> tau_iso = StandardTensors<dim>::dev_P * tau_bar;
  AssertThrow(std::fabs((symmetrize(F * static_cast<Tensor<2, dim>>(S_iso) *
                                    transpose(F)) -
                         tau_iso)
                          .norm()) < 1e-9,
              ExcInternalError());
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test_standard_tensors<2>();
  test_standard_tensors<3>();

  deallog << "OK" << std::endl;
}
