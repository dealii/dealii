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

#include <deal.II/base/polynomials_piecewise.h>


DEAL_II_NAMESPACE_OPEN



namespace Polynomials
{

  template <typename number>
  PiecewisePolynomial<number>::PiecewisePolynomial (const Polynomial<number> &coefficients_on_interval,
                                                    const unsigned int        n_intervals,
                                                    const unsigned int        interval,
                                                    const bool                spans_next_interval)
    :
    polynomial               (coefficients_on_interval),
    n_intervals              (n_intervals),
    interval                 (interval),
    spans_two_intervals      (spans_next_interval)
  {
    Assert (n_intervals > 0, ExcMessage ("No intervals given"));
    AssertIndexRange (interval, n_intervals);
  }



  template <typename number>
  void
  PiecewisePolynomial<number>::value (const number         x,
                                      std::vector<number> &values) const
  {
    Assert (values.size() > 0, ExcZero());
    const unsigned int values_size=values.size();

    // shift polynomial if necessary
    number y = x;
    double derivative_change_sign = 1.;
    if (n_intervals > 0)
      {
        const number step = 1./n_intervals;
        // polynomial spans over two intervals
        if (spans_two_intervals)
          {
            const double offset = step * interval;
            if (x<offset || x>offset+step+step)
              {
                for (unsigned int k=0; k<values.size(); ++k)
                  values[k] = 0;
                return;
              }
            else if (x<offset+step)
              y = x-offset;
            else
              {
                derivative_change_sign = -1.;
                y = offset+step+step-x;
              }
          }
        else
          {
            const double offset = step * interval;
            if (x<offset || x>offset+step)
              {
                for (unsigned int k=0; k<values.size(); ++k)
                  values[k] = 0;
                return;
              }
            else
              y = x-offset;
          }

        // on subinterval boundaries, cannot evaluate derivatives properly, so
        // set them to zero
        if ((std::abs(y)<1e-14 && (interval > 0 ||
                                   derivative_change_sign == -1.))
            ||
            (std::abs(y-step)<1e-14 &&
             (interval < n_intervals-1 || derivative_change_sign == -1.)))
          {
            values[0] = value(x);
            for (unsigned int d=1; d<values_size; ++d)
              values[d] = 0;
            return;
          }
      }

    polynomial.value(y, values);

    // change sign if necessary
    for (unsigned int j=1; j<values_size; j+=2)
      values[j] *= derivative_change_sign;
  }



  std::vector<PiecewisePolynomial<double> >
  generate_complete_Lagrange_basis_on_subdivisions (const unsigned int n_subdivisions,
                                                    const unsigned int base_degree)
  {
    std::vector<Polynomial<double> > p_base =
      LagrangeEquidistant::generate_complete_basis(base_degree);
    for (unsigned int i=0; i<p_base.size(); ++i)
      p_base[i].scale(n_subdivisions);

    std::vector<PiecewisePolynomial<double> > p;
    p.reserve (n_subdivisions * base_degree + 1);

    p.push_back (PiecewisePolynomial<double> (p_base[0], n_subdivisions, 0,
                                              false));
    for (unsigned int s=0; s<n_subdivisions; ++s)
      for (unsigned int i=0; i<base_degree; ++i)
        p.push_back (PiecewisePolynomial<double> (p_base[i+1], n_subdivisions,
                                                  s,
                                                  i==(base_degree-1) &&
                                                  s<n_subdivisions-1));
    return p;
  }

}

// ------------------ explicit instantiations --------------- //

namespace Polynomials
{
  template class PiecewisePolynomial<float>;
  template class PiecewisePolynomial<double>;
  template class PiecewisePolynomial<long double>;
}

DEAL_II_NAMESPACE_CLOSE
