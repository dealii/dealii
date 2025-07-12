// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Plot PolynomialsRT_Bubbles on the reference cell

#include <deal.II/base/polynomials_rt_bubbles.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <vector>

#include "../tests.h"

template <int dim>
std::vector<unsigned int>
get_lexicographic_numbering_RT(const unsigned int degree);

template <int dim>
void
plot(const PolynomialsRT_Bubbles<dim> &poly)
{
  QTrapezoid<1>               base_quadrature;
  QIterated<dim>              quadrature(base_quadrature, poly.degree() + 2);
  std::vector<Tensor<1, dim>> values(poly.n());
  std::vector<Tensor<2, dim>> grads;
  std::vector<Tensor<3, dim>> grads2;
  std::vector<Tensor<4, dim>> thirds;
  std::vector<Tensor<5, dim>> fourths;

  for (unsigned int k = 0; k < quadrature.size(); ++k)
    {
      if (k % (poly.degree() + 3) == 0)
        deallog << "RT_Bubbles" << poly.degree() << '<' << dim << '>'
                << std::endl;

      deallog << "RT_Bubbles" << poly.degree() << '<' << dim << '>' << '\t'
              << quadrature.point(k);
      poly.evaluate(
        quadrature.point(k), values, grads, grads2, thirds, fourths);

      for (unsigned int i = 0; i < poly.n(); ++i)
        for (unsigned int d = 0; d < dim; ++d)
          deallog << '\t' << values[i][d];
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  PolynomialsRT_Bubbles<2> p20(1, get_lexicographic_numbering_RT<2>(0));
  PolynomialsRT_Bubbles<2> p21(2, get_lexicographic_numbering_RT<2>(1));
  PolynomialsRT_Bubbles<2> p22(3, get_lexicographic_numbering_RT<2>(2));

  plot(p20);
  plot(p21);
  plot(p22);

  PolynomialsRT_Bubbles<3> p30(1, get_lexicographic_numbering_RT<3>(0));
  PolynomialsRT_Bubbles<3> p31(2, get_lexicographic_numbering_RT<3>(1));
  PolynomialsRT_Bubbles<3> p32(3, get_lexicographic_numbering_RT<3>(2));

  plot(p30);
  plot(p31);
  plot(p32);
}



template <int dim>
std::vector<unsigned int>
get_lexicographic_numbering_RT(const unsigned int degree)
{
  const unsigned int        n_dofs_face = Utilities::pow(degree + 1, dim - 1);
  std::vector<unsigned int> lexicographic_numbering;
  // component 1
  for (unsigned int j = 0; j < n_dofs_face; ++j)
    {
      lexicographic_numbering.push_back(j);
      if (degree + 1 > 1)
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

      if (degree + 1 > 1)
        for (unsigned int i = n_dofs_face * (2 * dim + degree);
             i < n_dofs_face * (2 * dim + degree) + degree * (degree + 1);
             ++i)
          {
            lexicographic_numbering.push_back(i + k_add * degree);
          }
      for (unsigned int j = n_dofs_face * 3; j < n_dofs_face * 3 + degree + 1;
           ++j)
        lexicographic_numbering.push_back(j + k_add);
    }

  // component 3
  if (dim == 3)
    {
      for (unsigned int i = 4 * n_dofs_face; i < 5 * n_dofs_face; ++i)
        lexicographic_numbering.push_back(i);
      if (degree + 1 > 1)
        for (unsigned int i = 6 * n_dofs_face + n_dofs_face * 2 * degree;
             i < 6 * n_dofs_face + n_dofs_face * 3 * degree;
             ++i)
          lexicographic_numbering.push_back(i);
      for (unsigned int i = 5 * n_dofs_face; i < 6 * n_dofs_face; ++i)
        lexicographic_numbering.push_back(i);
    }

  return lexicographic_numbering;
}
