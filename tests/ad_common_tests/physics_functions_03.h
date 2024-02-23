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
#include <deal.II/physics/transformations.h>

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
  const Tensor<1, dim, ADNumberType>          N{};
  const Tensor<1, dim, ADNumberType>          V{};
  const Tensor<2, dim, ADNumberType>          T1{};
  const SymmetricTensor<2, dim, ADNumberType> T2{};
  const Tensor<4, dim, ADNumberType>          TT1{};
  const SymmetricTensor<4, dim, ADNumberType> TT2{};

  const Tensor<1, dim, ADNumberType> ddet_F_dC =
    Physics::Transformations::nansons_formula(N, F);

  // Contravariant
  {
    const Tensor<1, dim, ADNumberType> V_out =
      Physics::Transformations::Contravariant::push_forward(V, F);
    const Tensor<2, dim, ADNumberType> T1_out =
      Physics::Transformations::Contravariant::push_forward(T1, F);
    const SymmetricTensor<2, dim, ADNumberType> T2_out =
      Physics::Transformations::Contravariant::push_forward(T2, F);
    const Tensor<4, dim, ADNumberType> TT1_out =
      Physics::Transformations::Contravariant::push_forward(TT1, F);
    const SymmetricTensor<4, dim, ADNumberType> TT2_out =
      Physics::Transformations::Contravariant::push_forward(TT2, F);
  }
  {
    const Tensor<1, dim, ADNumberType> V_out =
      Physics::Transformations::Contravariant::pull_back(V, F);
    const Tensor<2, dim, ADNumberType> T1_out =
      Physics::Transformations::Contravariant::pull_back(T1, F);
    const SymmetricTensor<2, dim, ADNumberType> T2_out =
      Physics::Transformations::Contravariant::pull_back(T2, F);
    const Tensor<4, dim, ADNumberType> TT1_out =
      Physics::Transformations::Contravariant::pull_back(TT1, F);
    const SymmetricTensor<4, dim, ADNumberType> TT2_out =
      Physics::Transformations::Contravariant::pull_back(TT2, F);
  }

  // Covariant
  {
    const Tensor<1, dim, ADNumberType> V_out =
      Physics::Transformations::Covariant::push_forward(V, F);
    const Tensor<2, dim, ADNumberType> T1_out =
      Physics::Transformations::Covariant::push_forward(T1, F);
    const SymmetricTensor<2, dim, ADNumberType> T2_out =
      Physics::Transformations::Covariant::push_forward(T2, F);
    const Tensor<4, dim, ADNumberType> TT1_out =
      Physics::Transformations::Covariant::push_forward(TT1, F);
    const SymmetricTensor<4, dim, ADNumberType> TT2_out =
      Physics::Transformations::Covariant::push_forward(TT2, F);
  }
  {
    const Tensor<1, dim, ADNumberType> V_out =
      Physics::Transformations::Covariant::pull_back(V, F);
    const Tensor<2, dim, ADNumberType> T1_out =
      Physics::Transformations::Covariant::pull_back(T1, F);
    const SymmetricTensor<2, dim, ADNumberType> T2_out =
      Physics::Transformations::Covariant::pull_back(T2, F);
    const Tensor<4, dim, ADNumberType> TT1_out =
      Physics::Transformations::Covariant::pull_back(TT1, F);
    const SymmetricTensor<4, dim, ADNumberType> TT2_out =
      Physics::Transformations::Covariant::pull_back(TT2, F);
  }

  // Piola
  {
    const Tensor<1, dim, ADNumberType> V_out =
      Physics::Transformations::Piola::push_forward(V, F);
    const Tensor<2, dim, ADNumberType> T1_out =
      Physics::Transformations::Piola::push_forward(T1, F);
    const SymmetricTensor<2, dim, ADNumberType> T2_out =
      Physics::Transformations::Piola::push_forward(T2, F);
    const Tensor<4, dim, ADNumberType> TT1_out =
      Physics::Transformations::Piola::push_forward(TT1, F);
    const SymmetricTensor<4, dim, ADNumberType> TT2_out =
      Physics::Transformations::Piola::push_forward(TT2, F);
  }
  {
    const Tensor<1, dim, ADNumberType> V_out =
      Physics::Transformations::Piola::pull_back(V, F);
    const Tensor<2, dim, ADNumberType> T1_out =
      Physics::Transformations::Piola::pull_back(T1, F);
    const SymmetricTensor<2, dim, ADNumberType> T2_out =
      Physics::Transformations::Piola::pull_back(T2, F);
    const Tensor<4, dim, ADNumberType> TT1_out =
      Physics::Transformations::Piola::pull_back(TT1, F);
    const SymmetricTensor<4, dim, ADNumberType> TT2_out =
      Physics::Transformations::Piola::pull_back(TT2, F);
  }
}
