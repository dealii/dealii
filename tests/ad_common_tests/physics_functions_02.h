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
// StandardTensors

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>

#include <iostream>

#include "../tests.h"

namespace AD = dealii::Differentiation::AD;

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_physics()
{
  using ADNumberType =
    typename AD::NumberTraits<number_t, ad_type_code>::ad_type;

  std::cout << "*** Test physics functions: Standard tensors, "
            << "dim = " << Utilities::to_string(dim) << ", "
            << "Type code: " << static_cast<int>(ad_type_code) << std::endl;

  Tensor<2, dim, ADNumberType>       grad_u;
  const Tensor<2, dim, ADNumberType> F =
    Physics::Elasticity::Kinematics::F(grad_u);

  const SymmetricTensor<2, dim, ADNumberType> ddet_F_dC =
    Physics::Elasticity::StandardTensors<dim>::ddet_F_dC(F);
  const SymmetricTensor<4, dim, ADNumberType> dC_inv_dC =
    Physics::Elasticity::StandardTensors<dim>::dC_inv_dC(F);

  const SymmetricTensor<4, dim, ADNumberType> Dev_P =
    Physics::Elasticity::StandardTensors<dim>::Dev_P(F);
  const SymmetricTensor<4, dim, ADNumberType> Dev_P_T =
    Physics::Elasticity::StandardTensors<dim>::Dev_P_T(F);
}
