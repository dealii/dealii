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

/**
 * plot Nedelec and Raviart-Thomas polynomials on the reference cell
 * to validate PolynomialsVectorAnisotropic.
 * We use the same numberings which are used by the actual elements.
 */

#include <deal.II/base/polynomials_vector_anisotropic.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <vector>

#include "../tests.h"

template <int dim>
std::vector<unsigned int>
get_lexicographic_numbering_RT(const unsigned int normal_degree,
                               const unsigned int tangential_degree);

template <int dim>
std::vector<unsigned int>
get_lexicographic_numbering_Nedelec(const unsigned int degree);

template <int dim>
void
plot(const PolynomialsVectorAnisotropic<dim> &poly)
{
  QTrapezoid<1>               base_quadrature;
  QIterated<dim>              quadrature(base_quadrature, poly.degree() + 3);
  std::vector<Tensor<1, dim>> values(poly.n());
  std::vector<Tensor<2, dim>> grads;
  std::vector<Tensor<3, dim>> grads2;
  std::vector<Tensor<4, dim>> thirds;
  std::vector<Tensor<5, dim>> fourths;



  for (unsigned int k = 0; k < quadrature.size(); ++k)
    {
      if (k % (poly.degree() + 4) == 0)
        deallog << "VectorAnisotropic" << poly.get_normal_degree()
                << poly.get_tangential_degree() << '<' << dim << '>'
                << std::endl;

      deallog << "VectorAnisotropic" << poly.get_normal_degree()
              << poly.get_tangential_degree() << '<' << dim << '>' << '\t'
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

  // Generating Nedelec polynomials in 2D
  PolynomialsVectorAnisotropic<2> p201(
    0, 1, get_lexicographic_numbering_Nedelec<2>(0));
  PolynomialsVectorAnisotropic<2> p212(
    1, 2, get_lexicographic_numbering_Nedelec<2>(1));
  PolynomialsVectorAnisotropic<2> p223(
    2, 3, get_lexicographic_numbering_Nedelec<2>(2));

  plot(p201);
  plot(p212);
  plot(p223);

  // Generating Nedelec polynomials in 3D
  PolynomialsVectorAnisotropic<3> p301(
    0, 1, get_lexicographic_numbering_Nedelec<3>(0));
  PolynomialsVectorAnisotropic<3> p312(
    1, 2, get_lexicographic_numbering_Nedelec<3>(1));
  PolynomialsVectorAnisotropic<3> p323(
    2, 3, get_lexicographic_numbering_Nedelec<3>(2));

  plot(p301);
  plot(p312);
  plot(p323);

  // Generating Raviart-Thomas polynomials in 2D
  PolynomialsVectorAnisotropic<2> p210(1,
                                       0,
                                       get_lexicographic_numbering_RT<2>(1, 0));
  PolynomialsVectorAnisotropic<2> p221(2,
                                       1,
                                       get_lexicographic_numbering_RT<2>(2, 1));
  PolynomialsVectorAnisotropic<2> p232(3,
                                       2,
                                       get_lexicographic_numbering_RT<2>(3, 2));

  plot(p210);
  plot(p221);
  plot(p232);

  // Generating Nedelec polynomials in 3D
  PolynomialsVectorAnisotropic<3> p310(1,
                                       0,
                                       get_lexicographic_numbering_RT<3>(1, 0));
  PolynomialsVectorAnisotropic<3> p321(2,
                                       1,
                                       get_lexicographic_numbering_RT<3>(2, 1));
  PolynomialsVectorAnisotropic<3> p332(3,
                                       2,
                                       get_lexicographic_numbering_RT<3>(3, 2));

  plot(p310);
  plot(p321);
  plot(p332);
}



template <int dim>
std::vector<unsigned int>
get_lexicographic_numbering_RT(const unsigned int normal_degree,
                               const unsigned int tangential_degree)
{
  const unsigned int n_dofs_face =
    Utilities::pow(tangential_degree + 1, dim - 1);
  std::vector<unsigned int> lexicographic_numbering;
  // component 1
  for (unsigned int j = 0; j < n_dofs_face; ++j)
    {
      lexicographic_numbering.push_back(j);
      if (normal_degree > 1)
        for (unsigned int i = n_dofs_face * 2 * dim;
             i < n_dofs_face * 2 * dim + normal_degree - 1;
             ++i)
          lexicographic_numbering.push_back(i + j * (normal_degree - 1));
      lexicographic_numbering.push_back(n_dofs_face + j);
    }

  // component 2
  unsigned int layers = (dim == 3) ? tangential_degree + 1 : 1;
  for (unsigned int k = 0; k < layers; ++k)
    {
      unsigned int k_add = k * (tangential_degree + 1);
      for (unsigned int j = n_dofs_face * 2;
           j < n_dofs_face * 2 + tangential_degree + 1;
           ++j)
        lexicographic_numbering.push_back(j + k_add);

      if (normal_degree > 1)
        for (unsigned int i = n_dofs_face * (2 * dim + (normal_degree - 1));
             i < n_dofs_face * (2 * dim + (normal_degree - 1)) +
                   (normal_degree - 1) * (tangential_degree + 1);
             ++i)
          {
            lexicographic_numbering.push_back(i + k_add * tangential_degree);
          }
      for (unsigned int j = n_dofs_face * 3;
           j < n_dofs_face * 3 + tangential_degree + 1;
           ++j)
        lexicographic_numbering.push_back(j + k_add);
    }

  // component 3
  if (dim == 3)
    {
      for (unsigned int i = 4 * n_dofs_face; i < 5 * n_dofs_face; ++i)
        lexicographic_numbering.push_back(i);
      if (normal_degree > 1)
        for (unsigned int i =
               6 * n_dofs_face + n_dofs_face * 2 * (normal_degree - 1);
             i < 6 * n_dofs_face + n_dofs_face * 3 * (normal_degree - 1);
             ++i)
          lexicographic_numbering.push_back(i);
      for (unsigned int i = 5 * n_dofs_face; i < 6 * n_dofs_face; ++i)
        lexicographic_numbering.push_back(i);
    }

  return lexicographic_numbering;
}



template <int dim>
std::vector<unsigned int>
get_lexicographic_numbering_Nedelec(const unsigned int degree)
{
  const unsigned int n_pols = (dim == 2) ?
                                (degree + 1) * (degree + 2) :
                                (degree + 1) * (degree + 2) * (degree + 2);

  std::vector<unsigned int> renumber_hierarchic_to_lexicographic;

  // edges 0-3
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back(n_pols + i * (degree + 2));
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back(n_pols + degree + 1 +
                                                   i * (degree + 2));
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back(i);
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back((degree + 1) * (degree + 1) +
                                                   i);
  if (dim == 2)
    {
      // quads for 2D
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(j * (degree + 1) + i);
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(n_pols +
                                                         j * (degree + 2) + i);
    }
  else if (dim == 3)
    {
      // edges 4-7
      const unsigned int offset_z = (degree + 1) * (degree + 1) * (degree + 2);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          n_pols + i * (degree + 2) + offset_z);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          n_pols + degree + 1 + i * (degree + 2) + offset_z);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(i + offset_z);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          (degree + 1) * (degree + 1) + i + offset_z);

      // edges 8-11
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + i * (degree + 2) * (degree + 2));
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + degree + 1 + i * (degree + 2) * (degree + 2));
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + (degree + 1) * (degree + 2) +
          i * (degree + 2) * (degree + 2));
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + (degree + 2) * (degree + 2) - 1 +
          i * (degree + 2) * (degree + 2));

      // quad 0
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            n_pols + i * (degree + 2) + j * (degree + 2) * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + i * (degree + 2) + j * (degree + 2) * (degree + 2));
      // quad 1
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            n_pols + degree + 1 + i * (degree + 2) +
            j * (degree + 2) * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + degree + 1 + i * (degree + 2) +
            j * (degree + 2) * (degree + 2));
      // quad 2
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(i + j * (degree + 2) *
                                                               (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + i + j * (degree + 2) * (degree + 2));
      // quad 3
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            (degree + 1) * (degree + 1) + i + j * (degree + 2) * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + (degree + 1) * (degree + 2) + i +
            j * (degree + 2) * (degree + 2));
      // quad 4
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(i + j * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(n_pols + i +
                                                         j * (degree + 2));
      // quad 5
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            (degree + 1) * (degree + 1) * (degree + 2) + i + j * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            n_pols + (degree + 1) * (degree + 1) * (degree + 2) + i +
            j * (degree + 2));

      // hexes
      for (unsigned int k = 1; k < degree + 1; ++k)
        for (unsigned int j = 1; j < degree + 1; ++j)
          for (unsigned int i = 0; i < degree + 1; ++i)
            renumber_hierarchic_to_lexicographic.push_back(
              k * (degree + 1) * (degree + 2) + j * (degree + 1) + i);
      for (unsigned int k = 1; k < degree + 1; ++k)
        for (unsigned int j = 0; j < degree + 1; ++j)
          for (unsigned int i = 1; i < degree + 1; ++i)
            renumber_hierarchic_to_lexicographic.push_back(
              n_pols + k * (degree + 2) * (degree + 1) + j * (degree + 2) + i);
      for (unsigned int k = 0; k < degree + 1; ++k)
        for (unsigned int j = 1; j < degree + 1; ++j)
          for (unsigned int i = 1; i < degree + 1; ++i)
            renumber_hierarchic_to_lexicographic.push_back(
              2 * n_pols + k * (degree + 2) * (degree + 2) + j * (degree + 2) +
              i);
    }

  return Utilities::invert_permutation(renumber_hierarchic_to_lexicographic);
}
