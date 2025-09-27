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

#ifndef dealii_polynomials_piecewise_h
#define dealii_polynomials_piecewise_h



#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Polynomials
 * @{
 */

/**
 * A namespace in which classes relating to the description of 1d polynomial
 * spaces are declared.
 */
namespace Polynomials
{
  /**
   * Definition of piecewise 1d polynomials for the unit interval. This space
   * allows the description of interpolating polynomials on parts of the unit
   * interval, similarly to the definition of finite element basis functions
   * on subdivided elements. The primary purpose of this class is to
   * allow constructing the shape functions of the FE_Q_iso_Q1 class that has
   * a number of interpolation points in each coordinate direction, but instead
   * of using them for higher-order polynomials just chooses piecewise linear
   * shape functions -- in effect, it is a $Q_1$ element defined on a
   * subdivision of the reference cell, and replicated on each of these
   * sub-cells.
   *
   * This class is not derived from the ScalarPolynomialsBase base class
   * because it is not actually a polynomial -- it is a piecewise polynomial.
   * However, it is interface-compatible with the Polynomials::Polynomial
   * class, and consequently can be used as template argument for
   * TensorProductPolynomials.
   *
   * @ingroup Polynomials
   */
  template <typename number>
  class PiecewisePolynomial : public EnableObserverPointer
  {
  public:
    /**
     * Constructor for Lagrange polynomial on an interval that is a subset of
     * the unit interval. It uses a polynomial description that is scaled to
     * the size of the subinterval compared to the unit interval, the total
     * number of intervals (subdivisions), the current index of the interval
     * as well as if the polynomial spans onto the next interval (e.g., if it
     * lives on two neighboring intervals).
     *
     * If the number of intervals is one, the piecewise polynomial behaves
     * exactly like a usual polynomial.
     */
    PiecewisePolynomial(const Polynomial<number> &coefficients_on_interval,
                        const unsigned int        n_intervals,
                        const unsigned int        interval,
                        const bool                spans_next_interval);

    /**
     * Constructor for linear Lagrange polynomial on an interval that is a
     * subset of the unit interval. It uses a polynomial description that is
     * scaled to the size of the subinterval compared to the unit interval.
     * The subintervals are bounded by the adjacent points in @p points.
     */
    PiecewisePolynomial(const std::vector<Point<1, number>> &points,
                        const unsigned int                   index);

    /**
     * Return the value of this polynomial at the given point, evaluating the
     * underlying polynomial. The polynomial evaluates to zero when outside of
     * the given interval (and possible the next one to the right when it
     * spans over that range).
     */
    number
    value(const number x) const;

    /**
     * Return the values and the derivatives of the Polynomial at point
     * <tt>x</tt>.  <tt>values[i], i=0,...,values.size()-1</tt> includes the
     * <tt>i</tt>th derivative. The number of derivatives to be computed is
     * thus determined by the size of the vector passed.
     *
     * Note that all the derivatives evaluate to zero at the border between
     * intervals (assuming exact arithmetic) in the interior of the unit
     * interval, as there is no unique gradient value in that case for a
     * piecewise polynomial. This is not always desired (e.g., when evaluating
     * jumps of gradients on the element boundary), but it is the user's
     * responsibility to avoid evaluation at these points when it does not
     * make sense.
     */
    void
    value(const number x, std::vector<number> &values) const;

    /**
     * Return the values and the derivatives of the Polynomial at point
     * <tt>x</tt>.  <tt>values[i], i=0,...,n_derivatives</tt> includes the
     * <tt>i</tt>th derivative.The number of derivatives to be computed is
     * determined by @p n_derivatives and @p values has to provide sufficient
     * space for @p n_derivatives + 1 values.
     *
     * Note that all the derivatives evaluate to zero at the border between
     * intervals (assuming exact arithmetic) in the interior of the unit
     * interval, as there is no unique gradient value in that case for a
     * piecewise polynomial. This is not always desired (e.g., when evaluating
     * jumps of gradients on the element boundary), but it is the user's
     * responsibility to avoid evaluation at these points when it does not
     * make sense.
     */
    void
    value(const number       x,
          const unsigned int n_derivatives,
          number            *values) const;

    /**
     * Degree of the polynomial. This is the degree of the underlying base
     * polynomial.
     */
    unsigned int
    degree() const;

    /**
     * Write or read the data of this object to or from a stream for the
     * purpose of serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version);

    /**
     * Return an estimate (in bytes) for the memory consumption of this object.
     */
    virtual std::size_t
    memory_consumption() const;

  protected:
    /**
     * Underlying polynomial object that is scaled to a subinterval and
     * concatenated accordingly.
     */
    Polynomial<number> polynomial;

    /**
     * A variable storing the number of intervals that the unit interval is
     * divided into.
     */
    unsigned int n_intervals;

    /**
     * A variable storing the index of the current polynomial in the range of
     * intervals.
     */
    unsigned int interval;

    /**
     * Store if the polynomial spans over two adjacent intervals, i.e., the
     * one given in subinterval and the next one.
     */
    bool spans_two_intervals;

    /**
     * Points bounding the subintervals in the case that piecewise linear
     * polynomial on varying subintervals is requested.
     */
    std::vector<number> points;

    /**
     * Precomputed inverses of the lengths of the subintervals, i.e.,
     * `one_over_lengths[i] = 1.0 / (points[i + 1] - points[i]`.
     */
    std::vector<number> one_over_lengths;

    /**
     * A variable storing the index of the current polynomial in the case that
     * piecewise linear polynomial on varying subintervals is requested.
     */
    unsigned int index;
  };



  /**
   * Generates a complete Lagrange basis on a subdivision of the unit interval
   * in smaller intervals for a given degree on the subintervals and number of
   * intervals.
   */
  std::vector<PiecewisePolynomial<double>>
  generate_complete_Lagrange_basis_on_subdivisions(
    const unsigned int n_subdivisions,
    const unsigned int base_degree);

  /**
   * Generates a complete linear basis on a subdivision of the unit interval
   * in smaller intervals for a given vector of points.
   */
  std::vector<PiecewisePolynomial<double>>
  generate_complete_linear_basis_on_subdivisions(
    const std::vector<Point<1>> &points);

} // namespace Polynomials


/** @} */

/* -------------------------- inline functions --------------------- */

namespace Polynomials
{
  template <typename number>
  inline unsigned int
  PiecewisePolynomial<number>::degree() const
  {
    if (points.size() > 0)
      return 1;
    return polynomial.degree();
  }



  template <typename number>
  inline number
  PiecewisePolynomial<number>::value(const number x) const
  {
    if (points.size() > 0)
      {
        if (x > points[index])
          return std::max<number>(0.0,
                                  1.0 - (x - points[index]) *
                                          one_over_lengths[index]);
        else if (x < points[index])
          return std::max<number>(0.0,
                                  0.0 + (x - points[index - 1]) *
                                          one_over_lengths[index - 1]);
        else
          return 1.0;
      }

    AssertIndexRange(interval, n_intervals);
    number y = x;
    // shift polynomial if necessary
    if (n_intervals > 1)
      {
        const number step = 1. / n_intervals;

        // polynomial spans over two intervals
        if (spans_two_intervals == true)
          {
            const number offset = step * interval;
            if (x < offset)
              return 0;
            else if (x > offset + step + step)
              return 0;
            else if (x < offset + step)
              y = x - offset;
            else
              y = offset + step + step - x;
          }
        else
          {
            const number offset = step * interval;
            if (x < offset || x > offset + step)
              return 0;
            else
              y = x - offset;
          }

        return polynomial.value(y);
      }
    else
      return polynomial.value(x);
  }



  template <typename number>
  template <class Archive>
  inline void
  PiecewisePolynomial<number>::serialize(Archive &ar, const unsigned int)
  {
    // forward to serialization function in the base class.
    ar &static_cast<EnableObserverPointer &>(*this);
    ar &polynomial;
    ar &n_intervals;
    ar &interval;
    ar &spans_two_intervals;
    ar &points;
    ar &one_over_lengths;
    ar &index;
  }
} // namespace Polynomials

DEAL_II_NAMESPACE_CLOSE

#endif
