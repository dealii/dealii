// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <int dim>
AutoDerivativeFunction<dim>::AutoDerivativeFunction(
  const double       hh,
  const unsigned int n_components,
  const double       initial_time)
  : Function<dim>(n_components, initial_time)
  , h(1)
  , ht(dim)
  , formula(Euler)
{
  set_h(hh);
  set_formula();
}



template <int dim>
void
AutoDerivativeFunction<dim>::set_formula(const DifferenceFormula form)
{
  // go through all known formulas, reject ones we don't know about
  // and don't handle in the member functions of this class
  switch (form)
    {
      case Euler:
      case UpwindEuler:
      case FourthOrder:
        break;
      default:
        Assert(false,
               ExcMessage("The argument passed to this function does not "
                          "match any known difference formula."));
    }

  formula = form;
}


template <int dim>
void
AutoDerivativeFunction<dim>::set_h(const double hh)
{
  h = hh;
  for (unsigned int i = 0; i < dim; ++i)
    ht[i][i] = h;
}


template <int dim>
Tensor<1, dim>
AutoDerivativeFunction<dim>::gradient(const Point<dim>  &p,
                                      const unsigned int comp) const
{
  Tensor<1, dim> grad;
  switch (formula)
    {
      case UpwindEuler:
        {
          Point<dim> q1;
          for (unsigned int i = 0; i < dim; ++i)
            {
              q1      = p - ht[i];
              grad[i] = (this->value(p, comp) - this->value(q1, comp)) / h;
            }
          break;
        }
      case Euler:
        {
          Point<dim> q1, q2;
          for (unsigned int i = 0; i < dim; ++i)
            {
              q1 = p + ht[i];
              q2 = p - ht[i];
              grad[i] =
                (this->value(q1, comp) - this->value(q2, comp)) / (2 * h);
            }
          break;
        }
      case FourthOrder:
        {
          Point<dim> q1, q2, q3, q4;
          for (unsigned int i = 0; i < dim; ++i)
            {
              q2      = p + ht[i];
              q1      = q2 + ht[i];
              q3      = p - ht[i];
              q4      = q3 - ht[i];
              grad[i] = (-this->value(q1, comp) + 8 * this->value(q2, comp) -
                         8 * this->value(q3, comp) + this->value(q4, comp)) /
                        (12 * h);
            }
          break;
        }
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
  return grad;
}


template <int dim>
void
AutoDerivativeFunction<dim>::vector_gradient(
  const Point<dim>            &p,
  std::vector<Tensor<1, dim>> &gradients) const
{
  Assert(gradients.size() == this->n_components,
         ExcDimensionMismatch(gradients.size(), this->n_components));

  switch (formula)
    {
      case UpwindEuler:
        {
          Point<dim>     q1;
          Vector<double> v(this->n_components), v1(this->n_components);
          const double   h_inv = 1. / h;
          for (unsigned int i = 0; i < dim; ++i)
            {
              q1 = p - ht[i];
              this->vector_value(p, v);
              this->vector_value(q1, v1);

              for (unsigned int comp = 0; comp < this->n_components; ++comp)
                gradients[comp][i] = (v(comp) - v1(comp)) * h_inv;
            }
          break;
        }

      case Euler:
        {
          Point<dim>     q1, q2;
          Vector<double> v1(this->n_components), v2(this->n_components);
          const double   h_inv_2 = 1. / (2 * h);
          for (unsigned int i = 0; i < dim; ++i)
            {
              q1 = p + ht[i];
              q2 = p - ht[i];
              this->vector_value(q1, v1);
              this->vector_value(q2, v2);

              for (unsigned int comp = 0; comp < this->n_components; ++comp)
                gradients[comp][i] = (v1(comp) - v2(comp)) * h_inv_2;
            }
          break;
        }

      case FourthOrder:
        {
          Point<dim>     q1, q2, q3, q4;
          Vector<double> v1(this->n_components), v2(this->n_components),
            v3(this->n_components), v4(this->n_components);
          const double h_inv_12 = 1. / (12 * h);
          for (unsigned int i = 0; i < dim; ++i)
            {
              q2 = p + ht[i];
              q1 = q2 + ht[i];
              q3 = p - ht[i];
              q4 = q3 - ht[i];
              this->vector_value(q1, v1);
              this->vector_value(q2, v2);
              this->vector_value(q3, v3);
              this->vector_value(q4, v4);

              for (unsigned int comp = 0; comp < this->n_components; ++comp)
                gradients[comp][i] =
                  (-v1(comp) + 8 * v2(comp) - 8 * v3(comp) + v4(comp)) *
                  h_inv_12;
            }
          break;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}


template <int dim>
void
AutoDerivativeFunction<dim>::gradient_list(
  const std::vector<Point<dim>> &points,
  std::vector<Tensor<1, dim>>   &gradients,
  const unsigned int             comp) const
{
  Assert(gradients.size() == points.size(),
         ExcDimensionMismatch(gradients.size(), points.size()));

  switch (formula)
    {
      case UpwindEuler:
        {
          Point<dim> q1;
          for (unsigned int p = 0; p < points.size(); ++p)
            for (unsigned int i = 0; i < dim; ++i)
              {
                q1 = points[p] - ht[i];
                gradients[p][i] =
                  (this->value(points[p], comp) - this->value(q1, comp)) / h;
              }
          break;
        }

      case Euler:
        {
          Point<dim> q1, q2;
          for (unsigned int p = 0; p < points.size(); ++p)
            for (unsigned int i = 0; i < dim; ++i)
              {
                q1 = points[p] + ht[i];
                q2 = points[p] - ht[i];
                gradients[p][i] =
                  (this->value(q1, comp) - this->value(q2, comp)) / (2 * h);
              }
          break;
        }

      case FourthOrder:
        {
          Point<dim> q1, q2, q3, q4;
          for (unsigned int p = 0; p < points.size(); ++p)
            for (unsigned int i = 0; i < dim; ++i)
              {
                q2 = points[p] + ht[i];
                q1 = q2 + ht[i];
                q3 = points[p] - ht[i];
                q4 = q3 - ht[i];
                gradients[p][i] =
                  (-this->value(q1, comp) + 8 * this->value(q2, comp) -
                   8 * this->value(q3, comp) + this->value(q4, comp)) /
                  (12 * h);
              }
          break;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim>
void
AutoDerivativeFunction<dim>::vector_gradient_list(
  const std::vector<Point<dim>>            &points,
  std::vector<std::vector<Tensor<1, dim>>> &gradients) const
{
  Assert(gradients.size() == points.size(),
         ExcDimensionMismatch(gradients.size(), points.size()));
  for (unsigned int p = 0; p < points.size(); ++p)
    Assert(gradients[p].size() == this->n_components,
           ExcDimensionMismatch(gradients.size(), this->n_components));

  switch (formula)
    {
      case UpwindEuler:
        {
          Point<dim> q1;
          for (unsigned int p = 0; p < points.size(); ++p)
            for (unsigned int i = 0; i < dim; ++i)
              {
                q1 = points[p] - ht[i];
                for (unsigned int comp = 0; comp < this->n_components; ++comp)
                  gradients[p][comp][i] =
                    (this->value(points[p], comp) - this->value(q1, comp)) / h;
              }
          break;
        }

      case Euler:
        {
          Point<dim> q1, q2;
          for (unsigned int p = 0; p < points.size(); ++p)
            for (unsigned int i = 0; i < dim; ++i)
              {
                q1 = points[p] + ht[i];
                q2 = points[p] - ht[i];
                for (unsigned int comp = 0; comp < this->n_components; ++comp)
                  gradients[p][comp][i] =
                    (this->value(q1, comp) - this->value(q2, comp)) / (2 * h);
              }
          break;
        }

      case FourthOrder:
        {
          Point<dim> q1, q2, q3, q4;
          for (unsigned int p = 0; p < points.size(); ++p)
            for (unsigned int i = 0; i < dim; ++i)
              {
                q2 = points[p] + ht[i];
                q1 = q2 + ht[i];
                q3 = points[p] - ht[i];
                q4 = q3 - ht[i];
                for (unsigned int comp = 0; comp < this->n_components; ++comp)
                  gradients[p][comp][i] =
                    (-this->value(q1, comp) + 8 * this->value(q2, comp) -
                     8 * this->value(q3, comp) + this->value(q4, comp)) /
                    (12 * h);
              }
          break;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}


template <int dim>
typename AutoDerivativeFunction<dim>::DifferenceFormula
AutoDerivativeFunction<dim>::get_formula_of_order(const unsigned int ord)
{
  switch (ord)
    {
      case 0:
      case 1:
        return UpwindEuler;
      case 2:
        return Euler;
      case 3:
      case 4:
        return FourthOrder;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
  return Euler;
}


template class AutoDerivativeFunction<1>;
template class AutoDerivativeFunction<2>;
template class AutoDerivativeFunction<3>;

DEAL_II_NAMESPACE_CLOSE
