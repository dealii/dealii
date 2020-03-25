// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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
  typedef
    typename AD::NumberTraits<number_t, ad_type_code>::ad_type ADNumberType;

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
