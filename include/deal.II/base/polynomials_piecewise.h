// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__polynomials_piecewise_h
#define __deal2__polynomials_piecewise_h



#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/point.h>

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
   * Definition of piecewise 1D polynomials for the unit interval. This space
   * allows the description of interpolating polynomials on parts of the unit
   * interval, similarly to the definition of finite element basis functions
   * on the subdivided elements. This primary purpose of this class is to
   * allow constructing FE_Q_iso_Q1 elements that put additional degrees of
   * freedom into an equivalent of a refined mesh instead of higher order
   * polynomials, which is useful when using mixed finite elements.
   *
   * @author Martin Kronbichler, 2013
   */
  template <typename number>
  class PiecewisePolynomial : public Subscriptor
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
    PiecewisePolynomial (const Polynomial<number> &coefficients_on_interval,
                         const unsigned int        n_intervals,
                         const unsigned int        interval,
                         const bool                spans_next_interval);

    /**
     * Return the value of this polynomial at the given point, evaluating the
     * underlying polynomial. The polynomial evaluates to zero when outside of
     * the given interval (and possible the next one to the right when it
     * spans over that range).
     */
    number value (const number x) const;

    /**
     * Return the values and the derivatives of the Polynomial at point
     * <tt>x</tt>.  <tt>values[i], i=0,...,values.size()-1</tt> includes the
     * <tt>i</tt>th derivative. The number of derivatives to be computed is
     * thus determined by the size of the array passed.
     *
     * Note that all the derivatives evaluate to zero at the border between
     * intervals (assuming exact arithmetics) in the interior of the unit
     * interval, as there is no unique gradient value in that case for a
     * piecewise polynomial. This is not always desired (e.g., when evaluating
     * jumps of gradients on the element boundary), but it is the user's
     * responsibility to avoid evaluation at these points when it does not
     * make sense.
     */
    void value (const number         x,
                std::vector<number> &values) const;

    /**
     * Degree of the polynomial. This is the degree of the underlying base
     * polynomial.
     */
    unsigned int degree () const;

    /**
     * Write or read the data of this object to or from a stream for the
     * purpose of serialization.
     */
    template <class Archive>
    void serialize (Archive &ar, const unsigned int version);

  protected:

    /**
     * Underlying polynomial object that is scaled to a subinterval and
     * concatenated accordingly.
     */
    Polynomial<number> polynomial;

    /**
     * Stores the number of intervals that the unit interval is divided into.
     */
    unsigned int n_intervals;

    /**
     * Stores the index of the current polynomial in the range of
     * intervals.
     */
    unsigned int interval;

    /**
     * Store if the polynomial spans over two adjacent intervals, i.e., the
     * one given in subinterval and the next one.
     */
    bool spans_two_intervals;
  };



  /**
   * Generates a complete Lagrange basis on a subdivision of the unit interval
   * in smaller intervals for a given degree on the subintervals and number of
   * intervals.
   */
  std::vector<PiecewisePolynomial<double> >
  generate_complete_Lagrange_basis_on_subdivisions (const unsigned int n_subdivisions,
                                                    const unsigned int base_degree);

}


/** @} */

/* -------------------------- inline functions --------------------- */

namespace Polynomials
{
  template <typename number>
  inline
  unsigned int
  PiecewisePolynomial<number>::degree () const
  {
    return polynomial.degree();
  }



  template <typename number>
  inline
  number
  PiecewisePolynomial<number>::value (const number x) const
  {
    AssertIndexRange (interval, n_intervals);
    number y = x;
    // shift polynomial if necessary
    if (n_intervals > 1)
      {
        const number step = 1./n_intervals;

        // polynomial spans over two intervals
        if (spans_two_intervals == true)
          {
            const number offset = step * interval;
            if (x<offset)
              return 0;
            else if (x>offset+step+step)
              return 0;
            else if (x<offset+step)
              y = x-offset;
            else
              y = offset+step+step-x;
          }
        else
          {
            const number offset = step * interval;
            if (x<offset || x>offset+step)
              return 0;
            else
              y = x-offset;
          }

        return polynomial.value(y);
      }
    else
      return polynomial.value(x);
  }



  template <typename number>
  template <class Archive>
  inline
  void
  PiecewisePolynomial<number>::serialize (Archive &ar, const unsigned int)
  {
    // forward to serialization function in the base class.
    ar &static_cast<Subscriptor &>(*this);
    ar &polynomial;
    ar &n_intervals;
    ar &interval;
    ar &spans_two_intervals;
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif
