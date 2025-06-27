// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/polynomials_raviart_thomas.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#include <iomanip>
#include <iostream>
#include <memory>


DEAL_II_NAMESPACE_OPEN



template <int dim>
PolynomialsRaviartThomas<dim>::PolynomialsRaviartThomas(const unsigned int k)
  : PolynomialsVectorAnisotropic<dim>(k + 1, k, get_lexicographic_numbering(k))
{}



template <int dim>
unsigned int
PolynomialsRaviartThomas<dim>::n_polynomials(const unsigned int degree)
{
  return PolynomialsVectorAnisotropic<dim>::n_polynomials(degree + 1, degree);
}



template <int dim>
std::vector<unsigned int>
PolynomialsRaviartThomas<dim>::get_lexicographic_numbering(
  const unsigned int degree)
{
  const unsigned int        n_dofs_face = Utilities::pow(degree + 1, dim - 1);
  std::vector<unsigned int> lexicographic_numbering;

  // component 1
  for (unsigned int j = 0; j < n_dofs_face; ++j)
    {
      lexicographic_numbering.push_back(j);
      for (unsigned int i = n_dofs_face * 2 * dim;
           i < n_dofs_face * 2 * dim + degree;
           ++i)
        lexicographic_numbering.push_back(i + j * degree);
      lexicographic_numbering.push_back(n_dofs_face + j);
    }

  // component 2
  unsigned int layers = (dim == 3) ? degree + 1 : 1;
  for (unsigned int k = 0; k < layers; ++k)
    {
      unsigned int k_add = k * (degree + 1);
      for (unsigned int j = n_dofs_face * 2; j < n_dofs_face * 2 + degree + 1;
           ++j)
        lexicographic_numbering.push_back(j + k_add);

      for (unsigned int i = n_dofs_face * (2 * dim + degree);
           i < n_dofs_face * (2 * dim + degree) + degree * (degree + 1);
           ++i)
        lexicographic_numbering.push_back(i + k_add * degree);
      for (unsigned int j = n_dofs_face * 3; j < n_dofs_face * 3 + degree + 1;
           ++j)
        lexicographic_numbering.push_back(j + k_add);
    }

  // component 3
  if (dim == 3)
    {
      for (unsigned int i = 4 * n_dofs_face; i < 5 * n_dofs_face; ++i)
        lexicographic_numbering.push_back(i);
      for (unsigned int i = 6 * n_dofs_face + n_dofs_face * 2 * degree;
           i < 6 * n_dofs_face + n_dofs_face * 3 * degree;
           ++i)
        lexicographic_numbering.push_back(i);
      for (unsigned int i = 5 * n_dofs_face; i < 6 * n_dofs_face; ++i)
        lexicographic_numbering.push_back(i);
    }

  return lexicographic_numbering;
}



template <int dim>
std::unique_ptr<TensorPolynomialsBase<dim>>
PolynomialsRaviartThomas<dim>::clone() const
{
  return std::make_unique<PolynomialsRaviartThomas<dim>>(*this);
}



template class PolynomialsRaviartThomas<1>;
template class PolynomialsRaviartThomas<2>;
template class PolynomialsRaviartThomas<3>;


DEAL_II_NAMESPACE_CLOSE
