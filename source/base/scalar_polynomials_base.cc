// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
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
#include <deal.II/base/scalar_polynomials_base.h>

#include <iomanip>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

template <int dim>
ScalarPolynomialsBase<dim>::ScalarPolynomialsBase(
  const unsigned int deg,
  const unsigned int n_polynomials)
  : polynomial_degree(deg)
  , n_pols(n_polynomials)
{
  // nothing to do here for now
}



template <int dim>
std::size_t
ScalarPolynomialsBase<dim>::memory_consumption() const
{
  DEAL_II_NOT_IMPLEMENTED();
  return 0;
}



template class ScalarPolynomialsBase<0>;
template class ScalarPolynomialsBase<1>;
template class ScalarPolynomialsBase<2>;
template class ScalarPolynomialsBase<3>;

DEAL_II_NAMESPACE_CLOSE
