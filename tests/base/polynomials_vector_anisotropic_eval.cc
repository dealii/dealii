// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/**
 * Check that PolynomialsVectorAnisotropic gives the same result irrespective
 * previous content in the result vector for values and derivatives.
 */

#include <deal.II/base/polynomials_vector_anisotropic.h>
#include <deal.II/base/tensor.h>

#include <vector>

#include "../tests.h"


template <int dim>
void
plot(const PolynomialsVectorAnisotropic<dim> &poly)
{
  // choose an evaluation point at some position within the reference cell
  Point<dim> eval_point;
  eval_point[0] = 0.371;
  if (dim > 1)
    eval_point[1] = 0.713;
  if (dim > 2)
    eval_point[2] = 0.137;

  std::vector<Tensor<1, dim>> values(poly.n());
  std::vector<Tensor<2, dim>> grads(poly.n());
  std::vector<Tensor<3, dim>> grads2(poly.n());
  std::vector<Tensor<4, dim>> thirds(poly.n());
  std::vector<Tensor<5, dim>> fourths(poly.n());


  deallog << "VectorAnisotropic" << poly.get_normal_degree()
          << poly.get_tangential_degree() << '<' << dim << '>' << " on point "
          << eval_point << std::endl;

  poly.evaluate(eval_point, values, grads, grads2, thirds, fourths);

  deallog << "values: " << std::endl;
  for (unsigned int i = 0; i < poly.n(); ++i)
    deallog << values[i] << std::endl;
  deallog << "grads: " << std::endl;
  for (unsigned int i = 0; i < poly.n(); ++i)
    deallog << grads[i] << std::endl;
  deallog << "grads2: " << std::endl;
  for (unsigned int i = 0; i < poly.n(); ++i)
    deallog << grads2[i] << std::endl;
  deallog << "thirds: " << std::endl;
  for (unsigned int i = 0; i < poly.n(); ++i)
    deallog << thirds[i] << std::endl;
  deallog << "fourths: " << std::endl;
  for (unsigned int i = 0; i < poly.n(); ++i)
    deallog << fourths[i] << std::endl;
  deallog << std::endl;

  deallog << "second evaluation" << std::endl;
  // fill some garbage content into the entries
  for (unsigned int i = 0; i < poly.n(); ++i)
    for (unsigned int d = 0; d < dim; ++d)
      values[i][d] = 0.2 * i + 0.33333 * d;
  for (unsigned int i = 0; i < poly.n(); ++i)
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int e = 0; e < dim; ++e)
        grads[i][d][e] = 0.2 * i + 0.33333 * d + e / 7.1;
  for (unsigned int i = 0; i < poly.n(); ++i)
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int e = 0; e < dim; ++e)
        for (unsigned int f = 0; f < dim; ++f)
          grads2[i][d][e][f] = 0.2 * i + 0.33333 * d + e / 7.1 - f * 19.2;
  for (unsigned int i = 0; i < poly.n(); ++i)
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int e = 0; e < dim; ++e)
        for (unsigned int f = 0; f < dim; ++f)
          for (unsigned int g = 0; g < dim; ++g)
            thirds[i][d][e][f][g] =
              0.2 * i + 0.33333 * d + e / 7.1 - f * 19.2 + g * 51.7;
  for (unsigned int i = 0; i < poly.n(); ++i)
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int e = 0; e < dim; ++e)
        for (unsigned int f = 0; f < dim; ++f)
          for (unsigned int g = 0; g < dim; ++g)
            for (unsigned int h = 0; h < dim; ++h)
              fourths[i][d][e][f][g][h] =
                0.2 * i + 0.33333 * d + e / 7.1 - f * 19.2 + g * 51.7 - 102 * h;
  poly.evaluate(eval_point, values, grads, grads2, thirds, fourths);

  deallog << "values: " << std::endl;
  for (unsigned int i = 0; i < poly.n(); ++i)
    deallog << values[i] << std::endl;
  deallog << "grads: " << std::endl;
  for (unsigned int i = 0; i < poly.n(); ++i)
    deallog << grads[i] << std::endl;
  deallog << "grads2: " << std::endl;
  for (unsigned int i = 0; i < poly.n(); ++i)
    deallog << grads2[i] << std::endl;
  deallog << "thirds: " << std::endl;
  for (unsigned int i = 0; i < poly.n(); ++i)
    deallog << thirds[i] << std::endl;
  deallog << "fourths: " << std::endl;
  for (unsigned int i = 0; i < poly.n(); ++i)
    deallog << fourths[i] << std::endl;
  deallog << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(8);

  // use Raviart-Thomas elements with some crazy numbering; use degree pair
  // (2,1) in 2D and (1,0) in 3D; the main polynomials are verified by other
  // tests, this test focuses on the second evaluation
  {
    const unsigned int        deg    = 1;
    const unsigned int        n_dofs = 2 * (deg + 1) * (deg + 2);
    std::vector<unsigned int> numbering(n_dofs);
    for (unsigned int i = 0; i < n_dofs / 2; ++i)
      numbering[i] = 2 * i;
    for (unsigned int i = 0; i < n_dofs / 2; ++i)
      numbering[n_dofs / 2 + i] = 2 * i + 1;
    PolynomialsVectorAnisotropic<2> pol(deg + 1, deg, numbering);
    plot(pol);
  }
  {
    const unsigned int        deg    = 0;
    const unsigned int        n_dofs = 3 * (deg + 1) * (deg + 1) * (deg + 2);
    std::vector<unsigned int> numbering(n_dofs);
    for (unsigned int i = 0; i < n_dofs / 3; ++i)
      numbering[i] = 3 * i;
    for (unsigned int i = 0; i < n_dofs / 3; ++i)
      numbering[n_dofs / 3 + i] = 3 * i + 1;
    for (unsigned int i = 0; i < n_dofs / 3; ++i)
      numbering[2 * n_dofs / 3 + i] = 3 * i + 2;
    PolynomialsVectorAnisotropic<3> pol(deg + 1, deg, numbering);
    plot(pol);
  }
}
