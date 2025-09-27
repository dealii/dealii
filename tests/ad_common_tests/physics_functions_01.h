// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Header file:
// Test to check that the functions in the Physics namespace compile
// with the various auto-differentiable number types:
// Kinematics

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/physics/elasticity/kinematics.h>

#include <iostream>

#include "../tests.h"

namespace AD = dealii::Differentiation::AD;

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_physics()
{
  using ADNumberType =
    typename AD::NumberTraits<number_t, ad_type_code>::ad_type;

  std::cout << "*** Test physics functions: Kinematics, "
            << "dim = " << Utilities::to_string(dim) << ", "
            << "Type code: " << static_cast<int>(ad_type_code) << std::endl;

  Tensor<2, dim, ADNumberType>       grad_u;
  const Tensor<2, dim, ADNumberType> F =
    Physics::Elasticity::Kinematics::F(grad_u);
  const Tensor<2, dim, ADNumberType> F_inv =
    invert(F); // Test, since this can be a problem too.

  const Tensor<2, dim, ADNumberType> F_iso =
    Physics::Elasticity::Kinematics::F_iso(F);
  const SymmetricTensor<2, dim, ADNumberType> F_vol =
    Physics::Elasticity::Kinematics::F_vol(F);
  const SymmetricTensor<2, dim, ADNumberType> C =
    Physics::Elasticity::Kinematics::C(F);
  const SymmetricTensor<2, dim, ADNumberType> b =
    Physics::Elasticity::Kinematics::b(F);

  const SymmetricTensor<2, dim, ADNumberType> E =
    Physics::Elasticity::Kinematics::E(F);
  const SymmetricTensor<2, dim, ADNumberType> epsilon =
    Physics::Elasticity::Kinematics::epsilon(grad_u);
  const SymmetricTensor<2, dim, ADNumberType> e =
    Physics::Elasticity::Kinematics::e(F);

  Tensor<2, dim, ADNumberType>       dF_dt;
  const Tensor<2, dim, ADNumberType> l =
    Physics::Elasticity::Kinematics::l(F, dF_dt);
  const SymmetricTensor<2, dim, ADNumberType> d =
    Physics::Elasticity::Kinematics::d(F, dF_dt);
  const Tensor<2, dim, ADNumberType> w =
    Physics::Elasticity::Kinematics::w(F, dF_dt);
}
