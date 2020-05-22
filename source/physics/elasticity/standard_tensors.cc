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
#include "standard_tensors.inst"

DEAL_II_NAMESPACE_CLOSE
