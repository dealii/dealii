// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
