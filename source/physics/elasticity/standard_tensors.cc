// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/physics/elasticity/standard_tensors.h>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN

template <int dim>
DEAL_II_CONSTEXPR const SymmetricTensor<2, dim>
                        Physics::Elasticity::StandardTensors<dim>::I
#  ifdef DEAL_II_CXX14_CONSTEXPR_BUG
  = unit_symmetric_tensor<dim>()
#  endif
  ;



template <int dim>
DEAL_II_CONSTEXPR const SymmetricTensor<4, dim>
                        Physics::Elasticity::StandardTensors<dim>::S
#  ifdef DEAL_II_CXX14_CONSTEXPR_BUG
  = identity_tensor<dim>()
#  endif
  ;



template <int dim>
DEAL_II_CONSTEXPR const SymmetricTensor<4, dim>
                        Physics::Elasticity::StandardTensors<dim>::IxI
#  ifdef DEAL_II_CXX14_CONSTEXPR_BUG
  = outer_product(unit_symmetric_tensor<dim>(), unit_symmetric_tensor<dim>())
#  endif
  ;



template <int dim>
DEAL_II_CONSTEXPR const SymmetricTensor<4, dim>
                        Physics::Elasticity::StandardTensors<dim>::dev_P
#  ifdef DEAL_II_CXX14_CONSTEXPR_BUG
  = deviator_tensor<dim>()
#  endif
  ;

#endif // DOXYGEN

// explicit instantiations
#include "physics/elasticity/standard_tensors.inst"

DEAL_II_NAMESPACE_CLOSE
