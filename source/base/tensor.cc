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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/config.h>

// Required for instantiation of template functions
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad/adolc_product_types.h>
#include <deal.II/differentiation/ad/sacado_product_types.h>

DEAL_II_NAMESPACE_OPEN


template <int dim, typename Number>
const unsigned int Tensor<0, dim, Number>::n_independent_components;

template <int rank, int dim, typename Number>
const unsigned int Tensor<rank, dim, Number>::n_independent_components;


#include "tensor.inst"


DEAL_II_NAMESPACE_CLOSE
