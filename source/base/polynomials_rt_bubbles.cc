// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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


#include "deal.II/base/polynomials_rt_bubbles.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#include <iomanip>
#include <iostream>
#include <memory>

DEAL_II_NAMESPACE_OPEN


template <int dim>
PolynomialsRT_Bubbles<dim>::PolynomialsRT_Bubbles(const unsigned int k)
  : TensorPolynomialsBase<dim>(k, n_polynomials(k))
  , raviart_thomas_space(k - 1)
  , monomials(k + 2)
{
  Assert(dim >= 2, ExcImpossibleInDim(dim));

  for (unsigned int i = 0; i < monomials.size(); ++i)
    monomials[i] = Polynomials::Monomial<double>(i);
}



template <int dim>
void
PolynomialsRT_Bubbles<dim>::evaluate(
  const Point<dim> &           unit_point,
  std::vector<Tensor<1, dim>> &values,
  std::vector<Tensor<2, dim>> &grads,
  std::vector<Tensor<3, dim>> &grad_grads,
  std::vector<Tensor<4, dim>> &third_derivatives,
  std::vector<Tensor<5, dim>> &fourth_derivatives) const
{
  Assert(values.size() == this->n() || values.size() == 0,
         ExcDimensionMismatch(values.size(), this->n()));
  Assert(grads.size() == this->n() || grads.size() == 0,
         ExcDimensionMismatch(grads.size(), this->n()));
  Assert(grad_grads.size() == this->n() || grad_grads.size() == 0,
         ExcDimensionMismatch(grad_grads.size(), this->n()));
  Assert(third_derivatives.size() == this->n() || third_derivatives.size() == 0,
         ExcDimensionMismatch(third_derivatives.size(), this->n()));
  Assert(fourth_derivatives.size() == this->n() ||
           fourth_derivatives.size() == 0,
         ExcDimensionMismatch(fourth_derivatives.size(), this->n()));

  // Third and fourth derivatives are not implemented
  (void)third_derivatives;
  Assert(third_derivatives.size() == 0, ExcNotImplemented());
  (void)fourth_derivatives;
  Assert(fourth_derivatives.size() == 0, ExcNotImplemented());

  const unsigned int n_sub     = raviart_thomas_space.n();
  const unsigned int my_degree = this->degree();

  // Guard access to the scratch arrays in the following block
  // using a mutex to make sure they are not used by multiple threads
  // at once
  {
    static std::mutex           mutex;
    std::lock_guard<std::mutex> lock(mutex);

    static std::vector<Tensor<1, dim>> p_values;
    static std::vector<Tensor<2, dim>> p_grads;
    static std::vector<Tensor<3, dim>> p_grad_grads;
    static std::vector<Tensor<4, dim>> p_third_derivatives;
    static std::vector<Tensor<5, dim>> p_fourth_derivatives;

    p_values.resize((values.size() == 0) ? 0 : n_sub);
    p_grads.resize((grads.size() == 0) ? 0 : n_sub);
    p_grad_grads.resize((grad_grads.size() == 0) ? 0 : n_sub);

    // This is the Raviart-Thomas part of the space
    raviart_thomas_space.evaluate(unit_point,
                                  p_values,
                                  p_grads,
                                  p_grad_grads,
                                  p_third_derivatives,
                                  p_fourth_derivatives);
    for (unsigned int i = 0; i < p_values.size(); ++i)
      values[i] = p_values[i];
    for (unsigned int i = 0; i < p_grads.size(); ++i)
      grads[i] = p_grads[i];
    for (unsigned int i = 0; i < p_grad_grads.size(); ++i)
      grad_grads[i] = p_grad_grads[i];
  }

  // Next we compute the polynomials and derivatives
  // of the curl part of the space
  const unsigned int n_derivatives = 3;
  double             monoval_plus[dim][n_derivatives + 1];
  double             monoval[dim][n_derivatives + 1];

  double monoval_i[dim][n_derivatives + 1];
  double monoval_j[dim][n_derivatives + 1];
  double monoval_jplus[dim][n_derivatives + 1];

  unsigned int start = n_sub;

  if (dim == 2)
    {
      // In 2d the curl part of the space is spanned by the vectors
      // of two types. The first one is
      //    [ x^i * [y^(k+1)]'  ]
      //    [ -[x^i]' * y^(k+1) ]
      // The second one can be obtained from the first by a cyclic
      // rotation of the coordinates.
      //  monoval_i = x^i,
      //  monoval_plus = x^(k+1)
      for (unsigned int d = 0; d < dim; ++d)
        monomials[my_degree + 1].value(unit_point(d),
                                       n_derivatives,
                                       monoval_plus[d]);

      for (unsigned int i = 0; i <= my_degree; ++i, ++start)
        {
          for (unsigned int d = 0; d < dim; ++d)
            monomials[i].value(unit_point(d), n_derivatives, monoval_i[d]);

          if (values.size() != 0)
            {
              values[start][0] = monoval_i[0][0] * monoval_plus[1][1];
              values[start][1] = -monoval_i[0][1] * monoval_plus[1][0];

              values[start + my_degree + 1][0] =
                -monoval_plus[0][0] * monoval_i[1][1];
              values[start + my_degree + 1][1] =
                monoval_plus[0][1] * monoval_i[1][0];
            }

          if (grads.size() != 0)
            {
              grads[start][0][0] = monoval_i[0][1] * monoval_plus[1][1];
              grads[start][0][1] = monoval_i[0][0] * monoval_plus[1][2];
              grads[start][1][0] = -monoval_i[0][2] * monoval_plus[1][0];
              grads[start][1][1] = -monoval_i[0][1] * monoval_plus[1][1];

              grads[start + my_degree + 1][0][0] =
                -monoval_plus[0][1] * monoval_i[1][1];
              grads[start + my_degree + 1][0][1] =
                -monoval_plus[0][0] * monoval_i[1][2];
              grads[start + my_degree + 1][1][0] =
                monoval_plus[0][2] * monoval_i[1][0];
              grads[start + my_degree + 1][1][1] =
                monoval_plus[0][1] * monoval_i[1][1];
            }

          if (grad_grads.size() != 0)
            {
              grad_grads[start][0][0][0] = monoval_i[0][2] * monoval_plus[1][1];
              grad_grads[start][0][0][1] = monoval_i[0][1] * monoval_plus[1][2];
              grad_grads[start][0][1][0] = monoval_i[0][1] * monoval_plus[1][2];
              grad_grads[start][0][1][1] = monoval_i[0][0] * monoval_plus[1][3];
              grad_grads[start][1][0][0] =
                -monoval_i[0][3] * monoval_plus[1][0];
              grad_grads[start][1][0][1] =
                -monoval_i[0][2] * monoval_plus[1][1];
              grad_grads[start][1][1][0] =
                -monoval_i[0][2] * monoval_plus[1][1];
              grad_grads[start][1][1][1] =
                -monoval_i[0][1] * monoval_plus[1][2];

              grad_grads[start + my_degree + 1][0][0][0] =
                -monoval_plus[0][2] * monoval_i[1][1];
              grad_grads[start + my_degree + 1][0][0][1] =
                -monoval_plus[0][1] * monoval_i[1][2];
              grad_grads[start + my_degree + 1][0][1][0] =
                -monoval_plus[0][1] * monoval_i[1][2];
              grad_grads[start + my_degree + 1][0][1][1] =
                -monoval_plus[0][0] * monoval_i[1][3];
              grad_grads[start + my_degree + 1][1][0][0] =
                monoval_plus[0][3] * monoval_i[1][0];
              grad_grads[start + my_degree + 1][1][0][1] =
                monoval_plus[0][2] * monoval_i[1][1];
              grad_grads[start + my_degree + 1][1][1][0] =
                monoval_plus[0][2] * monoval_i[1][1];
              grad_grads[start + my_degree + 1][1][1][1] =
                monoval_plus[0][1] * monoval_i[1][2];
            }
        }
      Assert(start == this->n() - my_degree - 1, ExcInternalError());
    }
  else if (dim == 3)
    {
      // In 3d the first type of basis vector is
      //  [ x^i * y^j * z^k * (j+k+2) ]
      //  [  -[x^i]' * y^(j+1) * z^k  ]
      //  [  -[x^i]' * y^j * z^(k+1)  ],
      // For the second type of basis vector y and z
      // are swapped. Then for each of these,
      // two more are obtained by the cyclic rotation
      // of the coordinates.
      //  monoval = x^k,   monoval_plus = x^(k+1)
      //  monoval_* = x^*, monoval_jplus = x^(j+1)
      for (unsigned int d = 0; d < dim; ++d)
        {
          monomials[my_degree + 1].value(unit_point(d),
                                         n_derivatives,
                                         monoval_plus[d]);
          monomials[my_degree].value(unit_point(d), n_derivatives, monoval[d]);
        }

      const unsigned int n_curls = (my_degree + 1) * (2 * my_degree + 1);
      // Span of $\tilde{B}$
      for (unsigned int i = 0; i <= my_degree; ++i)
        {
          for (unsigned int d = 0; d < dim; ++d)
            monomials[i].value(unit_point(d), n_derivatives, monoval_i[d]);

          for (unsigned int j = 0; j <= my_degree; ++j)
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  monomials[j].value(unit_point(d),
                                     n_derivatives,
                                     monoval_j[d]);
                  monomials[j + 1].value(unit_point(d),
                                         n_derivatives,
                                         monoval_jplus[d]);
                }

              if (values.size() != 0)
                {
                  values[start][0] = monoval_i[0][0] * monoval_j[1][0] *
                                     monoval[2][0] *
                                     static_cast<double>(j + my_degree + 2);
                  values[start][1] =
                    -monoval_i[0][1] * monoval_jplus[1][0] * monoval[2][0];
                  values[start][2] =
                    -monoval_i[0][1] * monoval_j[1][0] * monoval_plus[2][0];

                  values[start + n_curls][0] =
                    -monoval_jplus[0][0] * monoval_i[1][1] * monoval[2][0];
                  values[start + n_curls][1] =
                    monoval_j[0][0] * monoval_i[1][0] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  values[start + n_curls][2] =
                    -monoval_j[0][0] * monoval_i[1][1] * monoval_plus[2][0];

                  values[start + 2 * n_curls][0] =
                    -monoval_jplus[0][0] * monoval[1][0] * monoval_i[2][1];
                  values[start + 2 * n_curls][1] =
                    -monoval_j[0][0] * monoval_plus[1][0] * monoval_i[2][1];
                  values[start + 2 * n_curls][2] =
                    monoval_j[0][0] * monoval[1][0] * monoval_i[2][0] *
                    static_cast<double>(j + my_degree + 2);

                  // Only unique triples of powers (i j k)
                  // and (i k j) are allowed, 0 <= i,j <= k
                  if (j != my_degree)
                    {
                      values[start + 1][0] =
                        monoval_i[0][0] * monoval[1][0] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      values[start + 1][1] =
                        -monoval_i[0][1] * monoval_plus[1][0] * monoval_j[2][0];
                      values[start + 1][2] =
                        -monoval_i[0][1] * monoval[1][0] * monoval_jplus[2][0];

                      values[start + n_curls + 1][0] =
                        -monoval_plus[0][0] * monoval_i[1][1] * monoval_j[2][0];
                      values[start + n_curls + 1][1] =
                        monoval[0][0] * monoval_i[1][0] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      values[start + n_curls + 1][2] =
                        -monoval[0][0] * monoval_i[1][1] * monoval_jplus[2][0];

                      values[start + 2 * n_curls + 1][0] =
                        -monoval_plus[0][0] * monoval_j[1][0] * monoval_i[2][1];
                      values[start + 2 * n_curls + 1][1] =
                        -monoval[0][0] * monoval_jplus[1][0] * monoval_i[2][1];
                      values[start + 2 * n_curls + 1][2] =
                        monoval[0][0] * monoval_j[1][0] * monoval_i[2][0] *
                        static_cast<double>(j + my_degree + 2);
                    }
                }

              if (grads.size() != 0)
                {
                  grads[start][0][0] = monoval_i[0][1] * monoval_j[1][0] *
                                       monoval[2][0] *
                                       static_cast<double>(j + my_degree + 2);
                  grads[start][0][1] = monoval_i[0][0] * monoval_j[1][1] *
                                       monoval[2][0] *
                                       static_cast<double>(j + my_degree + 2);
                  grads[start][0][2] = monoval_i[0][0] * monoval_j[1][0] *
                                       monoval[2][1] *
                                       static_cast<double>(j + my_degree + 2);
                  grads[start][1][0] =
                    -monoval_i[0][2] * monoval_jplus[1][0] * monoval[2][0];
                  grads[start][1][1] =
                    -monoval_i[0][1] * monoval_jplus[1][1] * monoval[2][0];
                  grads[start][1][2] =
                    -monoval_i[0][1] * monoval_jplus[1][0] * monoval[2][1];
                  grads[start][2][0] =
                    -monoval_i[0][2] * monoval_j[1][0] * monoval_plus[2][0];
                  grads[start][2][1] =
                    -monoval_i[0][1] * monoval_j[1][1] * monoval_plus[2][0];
                  grads[start][2][2] =
                    -monoval_i[0][1] * monoval_j[1][0] * monoval_plus[2][1];

                  grads[start + n_curls][0][0] =
                    -monoval_jplus[0][1] * monoval_i[1][1] * monoval[2][0];
                  grads[start + n_curls][0][1] =
                    -monoval_jplus[0][0] * monoval_i[1][2] * monoval[2][0];
                  grads[start + n_curls][0][2] =
                    -monoval_jplus[0][0] * monoval_i[1][1] * monoval[2][1];
                  grads[start + n_curls][1][0] =
                    monoval_j[0][1] * monoval_i[1][0] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grads[start + n_curls][1][1] =
                    monoval_j[0][0] * monoval_i[1][1] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grads[start + n_curls][1][2] =
                    monoval_j[0][0] * monoval_i[1][0] * monoval[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grads[start + n_curls][2][0] =
                    -monoval_j[0][1] * monoval_i[1][1] * monoval_plus[2][0];
                  grads[start + n_curls][2][1] =
                    -monoval_j[0][0] * monoval_i[1][2] * monoval_plus[2][0];
                  grads[start + n_curls][2][2] =
                    -monoval_j[0][0] * monoval_i[1][1] * monoval_plus[2][1];

                  grads[start + 2 * n_curls][0][0] =
                    -monoval_jplus[0][1] * monoval[1][0] * monoval_i[2][1];
                  grads[start + 2 * n_curls][0][1] =
                    -monoval_jplus[0][0] * monoval[1][1] * monoval_i[2][1];
                  grads[start + 2 * n_curls][0][2] =
                    -monoval_jplus[0][0] * monoval[1][0] * monoval_i[2][2];
                  grads[start + 2 * n_curls][1][0] =
                    -monoval_j[0][1] * monoval_plus[1][0] * monoval_i[2][1];
                  grads[start + 2 * n_curls][1][1] =
                    -monoval_j[0][0] * monoval_plus[1][1] * monoval_i[2][1];
                  grads[start + 2 * n_curls][1][2] =
                    -monoval_j[0][0] * monoval_plus[1][0] * monoval_i[2][2];
                  grads[start + 2 * n_curls][2][0] =
                    monoval_j[0][1] * monoval[1][0] * monoval_i[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grads[start + 2 * n_curls][2][1] =
                    monoval_j[0][0] * monoval[1][1] * monoval_i[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grads[start + 2 * n_curls][2][2] =
                    monoval_j[0][0] * monoval[1][0] * monoval_i[2][1] *
                    static_cast<double>(j + my_degree + 2);

                  if (j != my_degree)
                    {
                      grads[start + 1][0][0] =
                        monoval_i[0][1] * monoval[1][0] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grads[start + 1][0][1] =
                        monoval_i[0][0] * monoval[1][1] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grads[start + 1][0][2] =
                        monoval_i[0][0] * monoval[1][0] * monoval_j[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grads[start + 1][1][0] =
                        -monoval_i[0][2] * monoval_plus[1][0] * monoval_j[2][0];
                      grads[start + 1][1][1] =
                        -monoval_i[0][1] * monoval_plus[1][1] * monoval_j[2][0];
                      grads[start + 1][1][2] =
                        -monoval_i[0][1] * monoval_plus[1][0] * monoval_j[2][1];
                      grads[start + 1][2][0] =
                        -monoval_i[0][2] * monoval[1][0] * monoval_jplus[2][0];
                      grads[start + 1][2][1] =
                        -monoval_i[0][1] * monoval[1][1] * monoval_jplus[2][0];
                      grads[start + 1][2][2] =
                        -monoval_i[0][1] * monoval[1][0] * monoval_jplus[2][1];

                      grads[start + n_curls + 1][0][0] =
                        -monoval_plus[0][1] * monoval_i[1][1] * monoval_j[2][0];
                      grads[start + n_curls + 1][0][1] =
                        -monoval_plus[0][0] * monoval_i[1][2] * monoval_j[2][0];
                      grads[start + n_curls + 1][0][2] =
                        -monoval_plus[0][0] * monoval_i[1][1] * monoval_j[2][1];
                      grads[start + n_curls + 1][1][0] =
                        monoval[0][1] * monoval_i[1][0] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grads[start + n_curls + 1][1][1] =
                        monoval[0][0] * monoval_i[1][1] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grads[start + n_curls + 1][1][2] =
                        monoval[0][0] * monoval_i[1][0] * monoval_j[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grads[start + n_curls + 1][2][0] =
                        -monoval[0][1] * monoval_i[1][1] * monoval_jplus[2][0];
                      grads[start + n_curls + 1][2][1] =
                        -monoval[0][0] * monoval_i[1][2] * monoval_jplus[2][0];
                      grads[start + n_curls + 1][2][2] =
                        -monoval[0][0] * monoval_i[1][1] * monoval_jplus[2][1];

                      grads[start + 2 * n_curls + 1][0][0] =
                        -monoval_plus[0][1] * monoval_j[1][0] * monoval_i[2][1];
                      grads[start + 2 * n_curls + 1][0][1] =
                        -monoval_plus[0][0] * monoval_j[1][1] * monoval_i[2][1];
                      grads[start + 2 * n_curls + 1][0][2] =
                        -monoval_plus[0][0] * monoval_j[1][0] * monoval_i[2][2];
                      grads[start + 2 * n_curls + 1][1][0] =
                        -monoval[0][1] * monoval_jplus[1][0] * monoval_i[2][1];
                      grads[start + 2 * n_curls + 1][1][1] =
                        -monoval[0][0] * monoval_jplus[1][1] * monoval_i[2][1];
                      grads[start + 2 * n_curls + 1][1][2] =
                        -monoval[0][0] * monoval_jplus[1][0] * monoval_i[2][2];
                      grads[start + 2 * n_curls + 1][2][0] =
                        monoval[0][1] * monoval_j[1][0] * monoval_i[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grads[start + 2 * n_curls + 1][2][1] =
                        monoval[0][0] * monoval_j[1][1] * monoval_i[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grads[start + 2 * n_curls + 1][2][2] =
                        monoval[0][0] * monoval_j[1][0] * monoval_i[2][1] *
                        static_cast<double>(j + my_degree + 2);
                    }
                }

              if (grad_grads.size() != 0)
                {
                  grad_grads[start][0][0][0] =
                    monoval_i[0][2] * monoval_j[1][0] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start][0][0][1] =
                    monoval_i[0][1] * monoval_j[1][1] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start][0][0][2] =
                    monoval_i[0][1] * monoval_j[1][0] * monoval[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start][0][1][0] =
                    monoval_i[0][1] * monoval_j[1][1] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start][0][1][1] =
                    monoval_i[0][0] * monoval_j[1][2] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start][0][1][2] =
                    monoval_i[0][0] * monoval_j[1][1] * monoval[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start][0][2][0] =
                    monoval_i[0][1] * monoval_j[1][0] * monoval[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start][0][2][1] =
                    monoval_i[0][0] * monoval_j[1][1] * monoval[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start][0][2][2] =
                    monoval_i[0][0] * monoval_j[1][0] * monoval[2][2] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start][1][0][0] =
                    -monoval_i[0][3] * monoval_jplus[1][0] * monoval[2][0];
                  grad_grads[start][1][0][1] =
                    -monoval_i[0][2] * monoval_jplus[1][1] * monoval[2][0];
                  grad_grads[start][1][0][2] =
                    -monoval_i[0][2] * monoval_jplus[1][0] * monoval[2][1];
                  grad_grads[start][1][1][0] =
                    -monoval_i[0][2] * monoval_jplus[1][1] * monoval[2][0];
                  grad_grads[start][1][1][1] =
                    -monoval_i[0][1] * monoval_jplus[1][2] * monoval[2][0];
                  grad_grads[start][1][1][2] =
                    -monoval_i[0][1] * monoval_jplus[1][1] * monoval[2][1];
                  grad_grads[start][1][2][0] =
                    -monoval_i[0][2] * monoval_jplus[1][0] * monoval[2][1];
                  grad_grads[start][1][2][1] =
                    -monoval_i[0][1] * monoval_jplus[1][1] * monoval[2][1];
                  grad_grads[start][1][2][2] =
                    -monoval_i[0][1] * monoval_jplus[1][0] * monoval[2][2];
                  grad_grads[start][2][0][0] =
                    -monoval_i[0][3] * monoval_j[1][0] * monoval_plus[2][0];
                  grad_grads[start][2][0][1] =
                    -monoval_i[0][2] * monoval_j[1][1] * monoval_plus[2][0];
                  grad_grads[start][2][0][2] =
                    -monoval_i[0][2] * monoval_j[1][0] * monoval_plus[2][1];
                  grad_grads[start][2][1][0] =
                    -monoval_i[0][2] * monoval_j[1][1] * monoval_plus[2][0];
                  grad_grads[start][2][1][1] =
                    -monoval_i[0][1] * monoval_j[1][2] * monoval_plus[2][0];
                  grad_grads[start][2][1][2] =
                    -monoval_i[0][1] * monoval_j[1][1] * monoval_plus[2][1];
                  grad_grads[start][2][2][0] =
                    -monoval_i[0][2] * monoval_j[1][0] * monoval_plus[2][1];
                  grad_grads[start][2][2][1] =
                    -monoval_i[0][1] * monoval_j[1][1] * monoval_plus[2][1];
                  grad_grads[start][2][2][2] =
                    -monoval_i[0][1] * monoval_j[1][0] * monoval_plus[2][2];

                  grad_grads[start + n_curls][0][0][0] =
                    -monoval_jplus[0][2] * monoval_i[1][1] * monoval[2][0];
                  grad_grads[start + n_curls][0][0][1] =
                    -monoval_jplus[0][1] * monoval_i[1][2] * monoval[2][0];
                  grad_grads[start + n_curls][0][0][2] =
                    -monoval_jplus[0][1] * monoval_i[1][1] * monoval[2][1];
                  grad_grads[start + n_curls][0][1][0] =
                    -monoval_jplus[0][1] * monoval_i[1][2] * monoval[2][0];
                  grad_grads[start + n_curls][0][1][1] =
                    -monoval_jplus[0][0] * monoval_i[1][3] * monoval[2][0];
                  grad_grads[start + n_curls][0][1][2] =
                    -monoval_jplus[0][0] * monoval_i[1][2] * monoval[2][1];
                  grad_grads[start + n_curls][0][2][0] =
                    -monoval_jplus[0][1] * monoval_i[1][1] * monoval[2][1];
                  grad_grads[start + n_curls][0][2][1] =
                    -monoval_jplus[0][0] * monoval_i[1][2] * monoval[2][1];
                  grad_grads[start + n_curls][0][2][2] =
                    -monoval_jplus[0][0] * monoval_i[1][1] * monoval[2][2];
                  grad_grads[start + n_curls][1][0][0] =
                    monoval_j[0][2] * monoval_i[1][0] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + n_curls][1][0][1] =
                    monoval_j[0][1] * monoval_i[1][1] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + n_curls][1][0][2] =
                    monoval_j[0][1] * monoval_i[1][0] * monoval[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + n_curls][1][1][0] =
                    monoval_j[0][1] * monoval_i[1][1] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + n_curls][1][1][1] =
                    monoval_j[0][0] * monoval_i[1][2] * monoval[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + n_curls][1][1][2] =
                    monoval_j[0][0] * monoval_i[1][1] * monoval[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + n_curls][1][2][0] =
                    monoval_j[0][1] * monoval_i[1][0] * monoval[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + n_curls][1][2][1] =
                    monoval_j[0][0] * monoval_i[1][1] * monoval[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + n_curls][1][2][2] =
                    monoval_j[0][0] * monoval_i[1][0] * monoval[2][2] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + n_curls][2][0][0] =
                    -monoval_j[0][2] * monoval_i[1][1] * monoval_plus[2][0];
                  grad_grads[start + n_curls][2][0][1] =
                    -monoval_j[0][1] * monoval_i[1][2] * monoval_plus[2][0];
                  grad_grads[start + n_curls][2][0][2] =
                    -monoval_j[0][1] * monoval_i[1][1] * monoval_plus[2][1];
                  grad_grads[start + n_curls][2][1][0] =
                    -monoval_j[0][1] * monoval_i[1][2] * monoval_plus[2][0];
                  grad_grads[start + n_curls][2][1][1] =
                    -monoval_j[0][0] * monoval_i[1][3] * monoval_plus[2][0];
                  grad_grads[start + n_curls][2][1][2] =
                    -monoval_j[0][0] * monoval_i[1][2] * monoval_plus[2][1];
                  grad_grads[start + n_curls][2][2][0] =
                    -monoval_j[0][1] * monoval_i[1][1] * monoval_plus[2][1];
                  grad_grads[start + n_curls][2][2][1] =
                    -monoval_j[0][0] * monoval_i[1][2] * monoval_plus[2][1];
                  grad_grads[start + n_curls][2][2][2] =
                    -monoval_j[0][0] * monoval_i[1][1] * monoval_plus[2][2];

                  grad_grads[start + 2 * n_curls][0][0][0] =
                    -monoval_jplus[0][2] * monoval[1][0] * monoval_i[2][1];
                  grad_grads[start + 2 * n_curls][0][0][1] =
                    -monoval_jplus[0][1] * monoval[1][1] * monoval_i[2][1];
                  grad_grads[start + 2 * n_curls][0][0][2] =
                    -monoval_jplus[0][1] * monoval[1][0] * monoval_i[2][2];
                  grad_grads[start + 2 * n_curls][0][1][0] =
                    -monoval_jplus[0][1] * monoval[1][1] * monoval_i[2][1];
                  grad_grads[start + 2 * n_curls][0][1][1] =
                    -monoval_jplus[0][0] * monoval[1][2] * monoval_i[2][1];
                  grad_grads[start + 2 * n_curls][0][1][2] =
                    -monoval_jplus[0][0] * monoval[1][1] * monoval_i[2][2];
                  grad_grads[start + 2 * n_curls][0][2][0] =
                    -monoval_jplus[0][1] * monoval[1][0] * monoval_i[2][2];
                  grad_grads[start + 2 * n_curls][0][2][1] =
                    -monoval_jplus[0][0] * monoval[1][1] * monoval_i[2][2];
                  grad_grads[start + 2 * n_curls][0][2][2] =
                    -monoval_jplus[0][0] * monoval[1][0] * monoval_i[2][3];
                  grad_grads[start + 2 * n_curls][1][0][0] =
                    -monoval_j[0][2] * monoval_plus[1][0] * monoval_i[2][1];
                  grad_grads[start + 2 * n_curls][1][0][1] =
                    -monoval_j[0][1] * monoval_plus[1][1] * monoval_i[2][1];
                  grad_grads[start + 2 * n_curls][1][0][2] =
                    -monoval_j[0][1] * monoval_plus[1][0] * monoval_i[2][2];
                  grad_grads[start + 2 * n_curls][1][1][0] =
                    -monoval_j[0][1] * monoval_plus[1][1] * monoval_i[2][1];
                  grad_grads[start + 2 * n_curls][1][1][1] =
                    -monoval_j[0][0] * monoval_plus[1][2] * monoval_i[2][1];
                  grad_grads[start + 2 * n_curls][1][1][2] =
                    -monoval_j[0][0] * monoval_plus[1][1] * monoval_i[2][2];
                  grad_grads[start + 2 * n_curls][1][2][0] =
                    -monoval_j[0][1] * monoval_plus[1][0] * monoval_i[2][2];
                  grad_grads[start + 2 * n_curls][1][2][1] =
                    -monoval_j[0][0] * monoval_plus[1][1] * monoval_i[2][2];
                  grad_grads[start + 2 * n_curls][1][2][2] =
                    -monoval_j[0][0] * monoval_plus[1][0] * monoval_i[2][3];
                  grad_grads[start + 2 * n_curls][2][0][0] =
                    monoval_j[0][2] * monoval[1][0] * monoval_i[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + 2 * n_curls][2][0][1] =
                    monoval_j[0][1] * monoval[1][1] * monoval_i[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + 2 * n_curls][2][0][2] =
                    monoval_j[0][1] * monoval[1][0] * monoval_i[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + 2 * n_curls][2][1][0] =
                    monoval_j[0][1] * monoval[1][1] * monoval_i[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + 2 * n_curls][2][1][1] =
                    monoval_j[0][0] * monoval[1][2] * monoval_i[2][0] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + 2 * n_curls][2][1][2] =
                    monoval_j[0][0] * monoval[1][1] * monoval_i[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + 2 * n_curls][2][2][0] =
                    monoval_j[0][1] * monoval[1][0] * monoval_i[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + 2 * n_curls][2][2][1] =
                    monoval_j[0][0] * monoval[1][1] * monoval_i[2][1] *
                    static_cast<double>(j + my_degree + 2);
                  grad_grads[start + 2 * n_curls][2][2][2] =
                    monoval_j[0][0] * monoval[1][0] * monoval_i[2][2] *
                    static_cast<double>(j + my_degree + 2);

                  if (j != my_degree)
                    {
                      grad_grads[start + 1][0][0][0] =
                        monoval_i[0][2] * monoval[1][0] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 1][0][0][1] =
                        monoval_i[0][1] * monoval[1][1] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 1][0][0][2] =
                        monoval_i[0][1] * monoval[1][0] * monoval_j[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 1][0][1][0] =
                        monoval_i[0][1] * monoval[1][1] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 1][0][1][1] =
                        monoval_i[0][0] * monoval[1][2] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 1][0][1][2] =
                        monoval_i[0][0] * monoval[1][1] * monoval_j[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 1][0][2][0] =
                        monoval_i[0][1] * monoval[1][0] * monoval_j[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 1][0][2][1] =
                        monoval_i[0][0] * monoval[1][1] * monoval_j[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 1][0][2][2] =
                        monoval_i[0][0] * monoval[1][0] * monoval_j[2][2] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 1][1][0][0] =
                        -monoval_i[0][3] * monoval_plus[1][0] * monoval_j[2][0];
                      grad_grads[start + 1][1][0][1] =
                        -monoval_i[0][2] * monoval_plus[1][1] * monoval_j[2][0];
                      grad_grads[start + 1][1][0][2] =
                        -monoval_i[0][2] * monoval_plus[1][0] * monoval_j[2][1];
                      grad_grads[start + 1][1][1][0] =
                        -monoval_i[0][2] * monoval_plus[1][1] * monoval_j[2][0];
                      grad_grads[start + 1][1][1][1] =
                        -monoval_i[0][1] * monoval_plus[1][2] * monoval_j[2][0];
                      grad_grads[start + 1][1][1][2] =
                        -monoval_i[0][1] * monoval_plus[1][1] * monoval_j[2][1];
                      grad_grads[start + 1][1][2][0] =
                        -monoval_i[0][2] * monoval_plus[1][0] * monoval_j[2][1];
                      grad_grads[start + 1][1][2][1] =
                        -monoval_i[0][1] * monoval_plus[1][1] * monoval_j[2][1];
                      grad_grads[start + 1][1][2][2] =
                        -monoval_i[0][1] * monoval_plus[1][0] * monoval_j[2][2];
                      grad_grads[start + 1][2][0][0] =
                        -monoval_i[0][3] * monoval[1][0] * monoval_jplus[2][0];
                      grad_grads[start + 1][2][0][1] =
                        -monoval_i[0][2] * monoval[1][1] * monoval_jplus[2][0];
                      grad_grads[start + 1][2][0][2] =
                        -monoval_i[0][2] * monoval[1][0] * monoval_jplus[2][1];
                      grad_grads[start + 1][2][1][0] =
                        -monoval_i[0][2] * monoval[1][1] * monoval_jplus[2][0];
                      grad_grads[start + 1][2][1][1] =
                        -monoval_i[0][1] * monoval[1][2] * monoval_jplus[2][0];
                      grad_grads[start + 1][2][1][2] =
                        -monoval_i[0][1] * monoval[1][1] * monoval_jplus[2][1];
                      grad_grads[start + 1][2][2][0] =
                        -monoval_i[0][2] * monoval[1][0] * monoval_jplus[2][1];
                      grad_grads[start + 1][2][2][1] =
                        -monoval_i[0][1] * monoval[1][1] * monoval_jplus[2][1];
                      grad_grads[start + 1][2][2][2] =
                        -monoval_i[0][1] * monoval[1][0] * monoval_jplus[2][2];

                      grad_grads[start + n_curls + 1][0][0][0] =
                        -monoval_plus[0][2] * monoval_i[1][1] * monoval_j[2][0];
                      grad_grads[start + n_curls + 1][0][0][1] =
                        -monoval_plus[0][1] * monoval_i[1][2] * monoval_j[2][0];
                      grad_grads[start + n_curls + 1][0][0][2] =
                        -monoval_plus[0][1] * monoval_i[1][1] * monoval_j[2][1];
                      grad_grads[start + n_curls + 1][0][1][0] =
                        -monoval_plus[0][1] * monoval_i[1][2] * monoval_j[2][0];
                      grad_grads[start + n_curls + 1][0][1][1] =
                        -monoval_plus[0][0] * monoval_i[1][3] * monoval_j[2][0];
                      grad_grads[start + n_curls + 1][0][1][2] =
                        -monoval_plus[0][0] * monoval_i[1][2] * monoval_j[2][1];
                      grad_grads[start + n_curls + 1][0][2][0] =
                        -monoval_plus[0][1] * monoval_i[1][1] * monoval_j[2][1];
                      grad_grads[start + n_curls + 1][0][2][1] =
                        -monoval_plus[0][0] * monoval_i[1][2] * monoval_j[2][1];
                      grad_grads[start + n_curls + 1][0][2][2] =
                        -monoval_plus[0][0] * monoval_i[1][1] * monoval_j[2][2];
                      grad_grads[start + n_curls + 1][1][0][0] =
                        monoval[0][2] * monoval_i[1][0] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + n_curls + 1][1][0][1] =
                        monoval[0][1] * monoval_i[1][1] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + n_curls + 1][1][0][2] =
                        monoval[0][1] * monoval_i[1][0] * monoval_j[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + n_curls + 1][1][1][0] =
                        monoval[0][1] * monoval_i[1][1] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + n_curls + 1][1][1][1] =
                        monoval[0][0] * monoval_i[1][2] * monoval_j[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + n_curls + 1][1][1][2] =
                        monoval[0][0] * monoval_i[1][1] * monoval_j[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + n_curls + 1][1][2][0] =
                        monoval[0][1] * monoval_i[1][0] * monoval_j[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + n_curls + 1][1][2][1] =
                        monoval[0][0] * monoval_i[1][1] * monoval_j[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + n_curls + 1][1][2][2] =
                        monoval[0][0] * monoval_i[1][0] * monoval_j[2][2] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + n_curls + 1][2][0][0] =
                        -monoval[0][2] * monoval_i[1][1] * monoval_jplus[2][0];
                      grad_grads[start + n_curls + 1][2][0][1] =
                        -monoval[0][1] * monoval_i[1][2] * monoval_jplus[2][0];
                      grad_grads[start + n_curls + 1][2][0][2] =
                        -monoval[0][1] * monoval_i[1][1] * monoval_jplus[2][1];
                      grad_grads[start + n_curls + 1][2][1][0] =
                        -monoval[0][1] * monoval_i[1][2] * monoval_jplus[2][0];
                      grad_grads[start + n_curls + 1][2][1][1] =
                        -monoval[0][0] * monoval_i[1][3] * monoval_jplus[2][0];
                      grad_grads[start + n_curls + 1][2][1][2] =
                        -monoval[0][0] * monoval_i[1][2] * monoval_jplus[2][1];
                      grad_grads[start + n_curls + 1][2][2][0] =
                        -monoval[0][1] * monoval_i[1][1] * monoval_jplus[2][1];
                      grad_grads[start + n_curls + 1][2][2][1] =
                        -monoval[0][0] * monoval_i[1][2] * monoval_jplus[2][1];
                      grad_grads[start + n_curls + 1][2][2][2] =
                        -monoval[0][0] * monoval_i[1][1] * monoval_jplus[2][2];

                      grad_grads[start + 2 * n_curls + 1][0][0][0] =
                        -monoval_plus[0][2] * monoval_j[1][0] * monoval_i[2][1];
                      grad_grads[start + 2 * n_curls + 1][0][0][1] =
                        -monoval_plus[0][1] * monoval_j[1][1] * monoval_i[2][1];
                      grad_grads[start + 2 * n_curls + 1][0][0][2] =
                        -monoval_plus[0][1] * monoval_j[1][0] * monoval_i[2][2];
                      grad_grads[start + 2 * n_curls + 1][0][1][0] =
                        -monoval_plus[0][1] * monoval_j[1][1] * monoval_i[2][1];
                      grad_grads[start + 2 * n_curls + 1][0][1][1] =
                        -monoval_plus[0][0] * monoval_j[1][2] * monoval_i[2][1];
                      grad_grads[start + 2 * n_curls + 1][0][1][2] =
                        -monoval_plus[0][0] * monoval_j[1][1] * monoval_i[2][2];
                      grad_grads[start + 2 * n_curls + 1][0][2][0] =
                        -monoval_plus[0][1] * monoval_j[1][0] * monoval_i[2][2];
                      grad_grads[start + 2 * n_curls + 1][0][2][1] =
                        -monoval_plus[0][0] * monoval_j[1][1] * monoval_i[2][2];
                      grad_grads[start + 2 * n_curls + 1][0][2][2] =
                        -monoval_plus[0][0] * monoval_j[1][0] * monoval_i[2][3];
                      grad_grads[start + 2 * n_curls + 1][1][0][0] =
                        -monoval[0][2] * monoval_jplus[1][0] * monoval_i[2][1];
                      grad_grads[start + 2 * n_curls + 1][1][0][1] =
                        -monoval[0][1] * monoval_jplus[1][1] * monoval_i[2][1];
                      grad_grads[start + 2 * n_curls + 1][1][0][2] =
                        -monoval[0][1] * monoval_jplus[1][0] * monoval_i[2][2];
                      grad_grads[start + 2 * n_curls + 1][1][1][0] =
                        -monoval[0][1] * monoval_jplus[1][1] * monoval_i[2][1];
                      grad_grads[start + 2 * n_curls + 1][1][1][1] =
                        -monoval[0][0] * monoval_jplus[1][2] * monoval_i[2][1];
                      grad_grads[start + 2 * n_curls + 1][1][1][2] =
                        -monoval[0][0] * monoval_jplus[1][1] * monoval_i[2][2];
                      grad_grads[start + 2 * n_curls + 1][1][2][0] =
                        -monoval[0][1] * monoval_jplus[1][0] * monoval_i[2][2];
                      grad_grads[start + 2 * n_curls + 1][1][2][1] =
                        -monoval[0][0] * monoval_jplus[1][1] * monoval_i[2][2];
                      grad_grads[start + 2 * n_curls + 1][1][2][2] =
                        -monoval[0][0] * monoval_jplus[1][0] * monoval_i[2][3];
                      grad_grads[start + 2 * n_curls + 1][2][0][0] =
                        monoval[0][2] * monoval_j[1][0] * monoval_i[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 2 * n_curls + 1][2][0][1] =
                        monoval[0][1] * monoval_j[1][1] * monoval_i[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 2 * n_curls + 1][2][0][2] =
                        monoval[0][1] * monoval_j[1][0] * monoval_i[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 2 * n_curls + 1][2][1][0] =
                        monoval[0][1] * monoval_j[1][1] * monoval_i[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 2 * n_curls + 1][2][1][1] =
                        monoval[0][0] * monoval_j[1][2] * monoval_i[2][0] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 2 * n_curls + 1][2][1][2] =
                        monoval[0][0] * monoval_j[1][1] * monoval_i[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 2 * n_curls + 1][2][2][0] =
                        monoval[0][1] * monoval_j[1][0] * monoval_i[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 2 * n_curls + 1][2][2][1] =
                        monoval[0][0] * monoval_j[1][1] * monoval_i[2][1] *
                        static_cast<double>(j + my_degree + 2);
                      grad_grads[start + 2 * n_curls + 1][2][2][2] =
                        monoval[0][0] * monoval_j[1][0] * monoval_i[2][2] *
                        static_cast<double>(j + my_degree + 2);
                    }
                }

              if (j == my_degree)
                start += 1;
              else
                start += 2;
            }
        }
      Assert(start == this->n() - 2 * n_curls, ExcInternalError());
    }
}



template <int dim>
unsigned int
PolynomialsRT_Bubbles<dim>::n_polynomials(const unsigned int k)
{
  if (dim == 1 || dim == 2 || dim == 3)
    return dim * Utilities::fixed_power<dim>(k + 1);

  Assert(false, ExcNotImplemented());
  return 0;
}


template <int dim>
std::unique_ptr<TensorPolynomialsBase<dim>>
PolynomialsRT_Bubbles<dim>::clone() const
{
  return std::make_unique<PolynomialsRT_Bubbles<dim>>(*this);
}


template class PolynomialsRT_Bubbles<1>;
template class PolynomialsRT_Bubbles<2>;
template class PolynomialsRT_Bubbles<3>;


DEAL_II_NAMESPACE_CLOSE
