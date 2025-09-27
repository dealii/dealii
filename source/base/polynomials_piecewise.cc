// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/polynomials_piecewise.h>


DEAL_II_NAMESPACE_OPEN



namespace Polynomials
{
  template <typename number>
  PiecewisePolynomial<number>::PiecewisePolynomial(
    const Polynomial<number> &coefficients_on_interval,
    const unsigned int        n_intervals,
    const unsigned int        interval,
    const bool                spans_next_interval)
    : polynomial(coefficients_on_interval)
    , n_intervals(n_intervals)
    , interval(interval)
    , spans_two_intervals(spans_next_interval)
    , index(numbers::invalid_unsigned_int)
  {
    Assert(n_intervals > 0, ExcMessage("No intervals given"));
    AssertIndexRange(interval, n_intervals);
  }



  template <typename number>
  PiecewisePolynomial<number>::PiecewisePolynomial(
    const std::vector<Point<1, number>> &points,
    const unsigned int                   index)
    : n_intervals(numbers::invalid_unsigned_int)
    , interval(numbers::invalid_unsigned_int)
    , spans_two_intervals(false)
    , index(index)
  {
    Assert(points.size() > 1, ExcMessage("No enough points given!"));
    AssertIndexRange(index, points.size());

    this->points.resize(points.size());
    for (unsigned int i = 0; i < points.size(); ++i)
      this->points[i] = points[i][0];

    this->one_over_lengths.resize(points.size() - 1);
    for (unsigned int i = 0; i < points.size() - 1; ++i)
      this->one_over_lengths[i] =
        number(1.0) / (points[i + 1][0] - points[i][0]);
  }



  template <typename number>
  void
  PiecewisePolynomial<number>::value(const number         x,
                                     std::vector<number> &values) const
  {
    Assert(values.size() > 0, ExcZero());

    value(x, values.size() - 1, values.data());
  }



  template <typename number>
  void
  PiecewisePolynomial<number>::value(const number       x,
                                     const unsigned int n_derivatives,
                                     number            *values) const
  {
    if (points.size() > 0)
      {
        if (x > points[index])
          values[0] = std::max<number>(0.0,
                                       1.0 - (x - points[index]) *
                                               one_over_lengths[index]);
        else if (x < points[index])
          values[0] = std::max<number>(0.0,
                                       0.0 + (x - points[index - 1]) *
                                               one_over_lengths[index - 1]);
        else
          values[0] = 1.0;

        if (n_derivatives >= 1)
          {
            if ((x > points[index]) && (points[index + 1] >= x))
              values[1] = -1.0 * one_over_lengths[index];
            else if ((x < points[index]) && (points[index - 1] <= x))
              values[1] = +1.0 * one_over_lengths[index - 1];
            else
              values[1] = 0.0;
          }

        // all other derivatives are zero
        for (unsigned int i = 2; i <= n_derivatives; ++i)
          values[i] = 0.0;

        return;
      }

    // shift polynomial if necessary
    number y                      = x;
    double derivative_change_sign = 1.;
    if (n_intervals > 0)
      {
        const number step = 1. / n_intervals;
        // polynomial spans over two intervals
        if (spans_two_intervals)
          {
            const double offset = step * interval;
            if (x < offset || x > offset + step + step)
              {
                for (unsigned int k = 0; k <= n_derivatives; ++k)
                  values[k] = 0;
                return;
              }
            else if (x < offset + step)
              y = x - offset;
            else
              {
                derivative_change_sign = -1.;
                y                      = offset + step + step - x;
              }
          }
        else
          {
            const double offset = step * interval;
            // ROCm 5.7 throws a floating point exception in debug mode when
            // trying to evaluate (x < offset || x > offset + step). Separating
            // the conditions fixes the issue.
            if (x < offset)
              {
                for (unsigned int k = 0; k <= n_derivatives; ++k)
                  values[k] = 0;
                return;
              }
            else if (x > offset + step)
              {
                for (unsigned int k = 0; k <= n_derivatives; ++k)
                  values[k] = 0;
                return;
              }
            else
              y = x - offset;
          }

        // on subinterval boundaries, cannot evaluate derivatives properly, so
        // set them to zero
        if ((std::abs(y) < 1e-14 &&
             (interval > 0 || derivative_change_sign == -1.)) ||
            (std::abs(y - step) < 1e-14 &&
             (interval < n_intervals - 1 || derivative_change_sign == -1.)))
          {
            values[0] = value(x);
            for (unsigned int d = 1; d <= n_derivatives; ++d)
              values[d] = 0;
            return;
          }
      }

    polynomial.value(y, n_derivatives, values);

    // change sign if necessary
    for (unsigned int j = 1; j <= n_derivatives; j += 2)
      values[j] *= derivative_change_sign;
  }



  template <typename number>
  std::size_t
  PiecewisePolynomial<number>::memory_consumption() const
  {
    return (polynomial.memory_consumption() +
            MemoryConsumption::memory_consumption(n_intervals) +
            MemoryConsumption::memory_consumption(interval) +
            MemoryConsumption::memory_consumption(spans_two_intervals) +
            MemoryConsumption::memory_consumption(points) +
            MemoryConsumption::memory_consumption(index));
  }



  std::vector<PiecewisePolynomial<double>>
  generate_complete_Lagrange_basis_on_subdivisions(
    const unsigned int n_subdivisions,
    const unsigned int base_degree)
  {
    std::vector<Polynomial<double>> p_base =
      LagrangeEquidistant::generate_complete_basis(base_degree);
    for (auto &polynomial : p_base)
      polynomial.scale(n_subdivisions);

    std::vector<PiecewisePolynomial<double>> p;
    p.reserve(n_subdivisions * base_degree + 1);

    p.emplace_back(p_base[0], n_subdivisions, 0, false);
    for (unsigned int s = 0; s < n_subdivisions; ++s)
      for (unsigned int i = 0; i < base_degree; ++i)
        p.emplace_back(p_base[i + 1],
                       n_subdivisions,
                       s,
                       i == (base_degree - 1) && s < n_subdivisions - 1);
    return p;
  }



  std::vector<PiecewisePolynomial<double>>
  generate_complete_linear_basis_on_subdivisions(
    const std::vector<Point<1>> &points)
  {
    std::vector<PiecewisePolynomial<double>> p;
    p.reserve(points.size());

    for (unsigned int s = 0; s < points.size(); ++s)
      p.emplace_back(points, s);

    return p;
  }

} // namespace Polynomials

// ------------------ explicit instantiations --------------- //

namespace Polynomials
{
  template class PiecewisePolynomial<float>;
  template class PiecewisePolynomial<double>;
  template class PiecewisePolynomial<long double>;
} // namespace Polynomials

DEAL_II_NAMESPACE_CLOSE
