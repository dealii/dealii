// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/table.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN


template <int dim>
unsigned int
PolynomialSpace<dim>::n_polynomials(const unsigned int n)
{
  unsigned int n_pols = n;
  for (unsigned int i = 1; i < dim; ++i)
    {
      n_pols *= (n + i);
      n_pols /= (i + 1);
    }
  return n_pols;
}


template <>
unsigned int
PolynomialSpace<0>::n_polynomials(const unsigned int)
{
  return 0;
}


template <>
std::array<unsigned int, 1>
PolynomialSpace<1>::compute_index(const unsigned int i) const
{
  AssertIndexRange(i, index_map.size());
  return {{index_map[i]}};
}



template <>
std::array<unsigned int, 2>
PolynomialSpace<2>::compute_index(const unsigned int i) const
{
  AssertIndexRange(i, index_map.size());
  const unsigned int n = index_map[i];
  // there should be a better way to
  // write this function (not
  // linear in n_1d), someone
  // should think about this...
  const unsigned int n_1d = polynomials.size();
  unsigned int       k    = 0;
  for (unsigned int iy = 0; iy < n_1d; ++iy)
    if (n < k + n_1d - iy)
      {
        return {{n - k, iy}};
      }
    else
      k += n_1d - iy;

  Assert(false, ExcInternalError());
  return {{numbers::invalid_unsigned_int, numbers::invalid_unsigned_int}};
}



template <>
std::array<unsigned int, 3>
PolynomialSpace<3>::compute_index(const unsigned int i) const
{
  AssertIndexRange(i, index_map.size());
  const unsigned int n = index_map[i];
  // there should be a better way to
  // write this function (not
  // quadratic in n_1d), someone
  // should think about this...
  //
  // (ah, and yes: the original
  // algorithm was even cubic!)
  const unsigned int n_1d = polynomials.size();
  unsigned int       k    = 0;
  for (unsigned int iz = 0; iz < n_1d; ++iz)
    for (unsigned int iy = 0; iy < n_1d - iz; ++iy)
      if (n < k + n_1d - iy - iz)
        {
          return {{n - k, iy, iz}};
        }
      else
        k += n_1d - iy - iz;

  Assert(false, ExcInternalError());
  return {{numbers::invalid_unsigned_int, numbers::invalid_unsigned_int}};
}


template <int dim>
void
PolynomialSpace<dim>::set_numbering(const std::vector<unsigned int> &renumber)
{
  Assert(renumber.size() == index_map.size(),
         ExcDimensionMismatch(renumber.size(), index_map.size()));

  index_map = renumber;
  for (unsigned int i = 0; i < index_map.size(); ++i)
    index_map_inverse[index_map[i]] = i;
}



template <int dim>
double
PolynomialSpace<dim>::compute_value(const unsigned int i,
                                    const Point<dim> & p) const
{
  const auto ix = compute_index(i);
  // take the product of the
  // polynomials in the various space
  // directions
  double result = 1.;
  for (unsigned int d = 0; d < dim; ++d)
    result *= polynomials[ix[d]].value(p(d));
  return result;
}



template <int dim>
Tensor<1, dim>
PolynomialSpace<dim>::compute_grad(const unsigned int i,
                                   const Point<dim> & p) const
{
  const auto ix = compute_index(i);

  Tensor<1, dim> result;
  for (unsigned int d = 0; d < dim; ++d)
    result[d] = 1.;

  // get value and first derivative
  std::vector<double> v(2);
  for (unsigned int d = 0; d < dim; ++d)
    {
      polynomials[ix[d]].value(p(d), v);
      result[d] *= v[1];
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        if (d1 != d)
          result[d1] *= v[0];
    }
  return result;
}


template <int dim>
Tensor<2, dim>
PolynomialSpace<dim>::compute_grad_grad(const unsigned int i,
                                        const Point<dim> & p) const
{
  const auto ix = compute_index(i);

  Tensor<2, dim> result;
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int d1 = 0; d1 < dim; ++d1)
      result[d][d1] = 1.;

  // get value, first and second
  // derivatives
  std::vector<double> v(3);
  for (unsigned int d = 0; d < dim; ++d)
    {
      polynomials[ix[d]].value(p(d), v);
      result[d][d] *= v[2];
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        {
          if (d1 != d)
            {
              result[d][d1] *= v[1];
              result[d1][d] *= v[1];
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                if (d2 != d)
                  result[d1][d2] *= v[0];
            }
        }
    }
  return result;
}


template <int dim>
void
PolynomialSpace<dim>::evaluate(
  const Point<dim> &           p,
  std::vector<double> &        values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  const unsigned int n_1d = polynomials.size();

  Assert(values.size() == this->n() || values.size() == 0,
         ExcDimensionMismatch2(values.size(), this->n(), 0));
  Assert(grads.size() == this->n() || grads.size() == 0,
         ExcDimensionMismatch2(grads.size(), this->n(), 0));
  Assert(grad_grads.size() == this->n() || grad_grads.size() == 0,
         ExcDimensionMismatch2(grad_grads.size(), this->n(), 0));
  Assert(third_derivatives.size() == this->n() || third_derivatives.size() == 0,
         ExcDimensionMismatch2(third_derivatives.size(), this->n(), 0));
  Assert(fourth_derivatives.size() == this->n() ||
           fourth_derivatives.size() == 0,
         ExcDimensionMismatch2(fourth_derivatives.size(), this->n(), 0));

  unsigned int v_size = 0;
  bool update_values = false, update_grads = false, update_grad_grads = false;
  bool update_3rd_derivatives = false, update_4th_derivatives = false;
  if (values.size() == this->n())
    {
      update_values = true;
      v_size        = 1;
    }
  if (grads.size() == this->n())
    {
      update_grads = true;
      v_size       = 2;
    }
  if (grad_grads.size() == this->n())
    {
      update_grad_grads = true;
      v_size            = 3;
    }
  if (third_derivatives.size() == this->n())
    {
      update_3rd_derivatives = true;
      v_size                 = 4;
    }
  if (fourth_derivatives.size() == this->n())
    {
      update_4th_derivatives = true;
      v_size                 = 5;
    }

  // Store data in a single
  // object. Access is by
  // v[d][n][o]
  //  d: coordinate direction
  //  n: number of 1d polynomial
  //  o: order of derivative
  Table<2, std::vector<double>> v(dim, n_1d);
  for (unsigned int d = 0; d < v.size()[0]; ++d)
    for (unsigned int i = 0; i < v.size()[1]; ++i)
      {
        v(d, i).resize(v_size, 0.);
        polynomials[i].value(p(d), v(d, i));
      }

  if (update_values)
    {
      unsigned int k = 0;

      for (unsigned int iz = 0; iz < ((dim > 2) ? n_1d : 1); ++iz)
        for (unsigned int iy = 0; iy < ((dim > 1) ? n_1d - iz : 1); ++iy)
          for (unsigned int ix = 0; ix < n_1d - iy - iz; ++ix)
            values[index_map_inverse[k++]] = v[0][ix][0] *
                                             ((dim > 1) ? v[1][iy][0] : 1.) *
                                             ((dim > 2) ? v[2][iz][0] : 1.);
    }

  if (update_grads)
    {
      unsigned int k = 0;

      for (unsigned int iz = 0; iz < ((dim > 2) ? n_1d : 1); ++iz)
        for (unsigned int iy = 0; iy < ((dim > 1) ? n_1d - iz : 1); ++iy)
          for (unsigned int ix = 0; ix < n_1d - iy - iz; ++ix)
            {
              const unsigned int k2 = index_map_inverse[k++];
              for (unsigned int d = 0; d < dim; ++d)
                grads[k2][d] = v[0][ix][(d == 0) ? 1 : 0] *
                               ((dim > 1) ? v[1][iy][(d == 1) ? 1 : 0] : 1.) *
                               ((dim > 2) ? v[2][iz][(d == 2) ? 1 : 0] : 1.);
            }
    }

  if (update_grad_grads)
    {
      unsigned int k = 0;

      for (unsigned int iz = 0; iz < ((dim > 2) ? n_1d : 1); ++iz)
        for (unsigned int iy = 0; iy < ((dim > 1) ? n_1d - iz : 1); ++iy)
          for (unsigned int ix = 0; ix < n_1d - iy - iz; ++ix)
            {
              const unsigned int k2 = index_map_inverse[k++];
              for (unsigned int d1 = 0; d1 < dim; ++d1)
                for (unsigned int d2 = 0; d2 < dim; ++d2)
                  {
                    // Derivative
                    // order for each
                    // direction
                    const unsigned int j0 =
                      ((d1 == 0) ? 1 : 0) + ((d2 == 0) ? 1 : 0);
                    const unsigned int j1 =
                      ((d1 == 1) ? 1 : 0) + ((d2 == 1) ? 1 : 0);
                    const unsigned int j2 =
                      ((d1 == 2) ? 1 : 0) + ((d2 == 2) ? 1 : 0);

                    grad_grads[k2][d1][d2] = v[0][ix][j0] *
                                             ((dim > 1) ? v[1][iy][j1] : 1.) *
                                             ((dim > 2) ? v[2][iz][j2] : 1.);
                  }
            }
    }

  if (update_3rd_derivatives)
    {
      unsigned int k = 0;

      for (unsigned int iz = 0; iz < ((dim > 2) ? n_1d : 1); ++iz)
        for (unsigned int iy = 0; iy < ((dim > 1) ? n_1d - iz : 1); ++iy)
          for (unsigned int ix = 0; ix < n_1d - iy - iz; ++ix)
            {
              const unsigned int k2 = index_map_inverse[k++];
              for (unsigned int d1 = 0; d1 < dim; ++d1)
                for (unsigned int d2 = 0; d2 < dim; ++d2)
                  for (unsigned int d3 = 0; d3 < dim; ++d3)
                    {
                      // Derivative
                      // order for each
                      // direction
                      std::vector<unsigned int> deriv_order(dim, 0);
                      for (unsigned int x = 0; x < dim; ++x)
                        {
                          if (d1 == x)
                            ++deriv_order[x];
                          if (d2 == x)
                            ++deriv_order[x];
                          if (d3 == x)
                            ++deriv_order[x];
                        }

                      third_derivatives[k2][d1][d2][d3] =
                        v[0][ix][deriv_order[0]] *
                        ((dim > 1) ? v[1][iy][deriv_order[1]] : 1.) *
                        ((dim > 2) ? v[2][iz][deriv_order[2]] : 1.);
                    }
            }
    }

  if (update_4th_derivatives)
    {
      unsigned int k = 0;

      for (unsigned int iz = 0; iz < ((dim > 2) ? n_1d : 1); ++iz)
        for (unsigned int iy = 0; iy < ((dim > 1) ? n_1d - iz : 1); ++iy)
          for (unsigned int ix = 0; ix < n_1d - iy - iz; ++ix)
            {
              const unsigned int k2 = index_map_inverse[k++];
              for (unsigned int d1 = 0; d1 < dim; ++d1)
                for (unsigned int d2 = 0; d2 < dim; ++d2)
                  for (unsigned int d3 = 0; d3 < dim; ++d3)
                    for (unsigned int d4 = 0; d4 < dim; ++d4)
                      {
                        // Derivative
                        // order for each
                        // direction
                        std::vector<unsigned int> deriv_order(dim, 0);
                        for (unsigned int x = 0; x < dim; ++x)
                          {
                            if (d1 == x)
                              ++deriv_order[x];
                            if (d2 == x)
                              ++deriv_order[x];
                            if (d3 == x)
                              ++deriv_order[x];
                            if (d4 == x)
                              ++deriv_order[x];
                          }

                        fourth_derivatives[k2][d1][d2][d3][d4] =
                          v[0][ix][deriv_order[0]] *
                          ((dim > 1) ? v[1][iy][deriv_order[1]] : 1.) *
                          ((dim > 2) ? v[2][iz][deriv_order[2]] : 1.);
                      }
            }
    }
}



template <int dim>
std::unique_ptr<ScalarPolynomialsBase<dim>>
PolynomialSpace<dim>::clone() const
{
  return std::make_unique<PolynomialSpace<dim>>(*this);
}


template class PolynomialSpace<1>;
template class PolynomialSpace<2>;
template class PolynomialSpace<3>;

DEAL_II_NAMESPACE_CLOSE
