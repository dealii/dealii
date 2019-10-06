// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_polynomials_base.h>

#include <iomanip>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

template <int dim>
TensorPolynomialsBase<dim>::TensorPolynomialsBase(
  const unsigned int deg,
  const unsigned int n_polynomials)
  : polynomial_degree(deg)
  , n_pols(n_polynomials)
{
  // nothing to do here for now
}



template class TensorPolynomialsBase<1>;
template class TensorPolynomialsBase<2>;
template class TensorPolynomialsBase<3>;

DEAL_II_NAMESPACE_CLOSE
