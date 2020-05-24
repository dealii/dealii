// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

DEAL_II_NAMESPACE_OPEN



// have a lock that guarantees that at most one thread is changing and
// accessing the @p{coefficients} arrays of classes implementing
// polynomials with tables. make this lock local to this file.
//
// having only one lock for all of these classes is probably not going
// to be a problem since we only need it on very rare occasions. if
// someone finds this is a bottleneck, feel free to replace it by a
// more fine-grained solution
namespace
{
  std::mutex coefficients_lock;
}



namespace Polynomials
{
  // -------------------- class Polynomial ---------------- //


  template <typename number>
  Polynomial<number>::Polynomial(const std::vector<number> &a)
    : coefficients(a)
    , in_lagrange_product_form(false)
    , lagrange_weight(1.)
  {}



  template <typename number>
  Polynomial<number>::Polynomial(const unsigned int n)
    : coefficients(n + 1, 0.)
    , in_lagrange_product_form(false)
    , lagrange_weight(1.)
  {}



  template <typename number>
  Polynomial<number>::Polynomial(const std::vector<Point<1>> &supp,
                                 const unsigned int           center)
    : in_lagrange_product_form(true)
  {
    Assert(supp.size() > 0, ExcEmptyObject());
    AssertIndexRange(center, supp.size());

    lagrange_support_points.reserve(supp.size() - 1);
    number tmp_lagrange_weight = 1.;
    for (unsigned int i = 0; i < supp.size(); ++i)
      if (i != center)
        {
          lagrange_support_points.push_back(supp[i](0));
          tmp_lagrange_weight *= supp[center](0) - supp[i](0);
        }

    // check for underflow and overflow
    Assert(std::fabs(tmp_lagrange_weight) > std::numeric_limits<number>::min(),
           ExcMessage("Underflow in computation of Lagrange denominator."));
    Assert(std::fabs(tmp_lagrange_weight) < std::numeric_limits<number>::max(),
           ExcMessage("Overflow in computation of Lagrange denominator."));

    lagrange_weight = static_cast<number>(1.) / tmp_lagrange_weight;
  }



  template <typename number>
  void
  Polynomial<number>::value(const number x, std::vector<number> &values) const
  {
    Assert(values.size() > 0, ExcZero());

    value(x, values.size() - 1, values.data());
  }



  template <typename number>
  void
  Polynomial<number>::value(const number       x,
                            const unsigned int n_derivatives,
                            number *           values) const
  {
    // evaluate Lagrange polynomial and derivatives
    if (in_lagrange_product_form == true)
      {
        // to compute the value and all derivatives of a polynomial of the
        // form (x-x_1)*(x-x_2)*...*(x-x_n), expand the derivatives like
        // automatic differentiation does.
        const unsigned int n_supp = lagrange_support_points.size();
        switch (n_derivatives)
          {
            default:
              values[0] = 1;
              for (unsigned int d = 1; d <= n_derivatives; ++d)
                values[d] = 0;
              for (unsigned int i = 0; i < n_supp; ++i)
                {
                  const number v = x - lagrange_support_points[i];

                  // multiply by (x-x_i) and compute action on all derivatives,
                  // too (inspired from automatic differentiation: implement the
                  // product rule for the old value and the new variable 'v',
                  // i.e., expand value v and derivative one). since we reuse a
                  // value from the next lower derivative from the steps before,
                  // need to start from the highest derivative
                  for (unsigned int k = n_derivatives; k > 0; --k)
                    values[k] = (values[k] * v + values[k - 1]);
                  values[0] *= v;
                }
              // finally, multiply by the weight in the Lagrange
              // denominator. Could be done instead of setting values[0] = 1
              // above, but that gives different accumulation of round-off
              // errors (multiplication is not associative) compared to when we
              // computed the weight, and hence a basis function might not be
              // exactly one at the center point, which is nice to have. We also
              // multiply derivatives by k! to transform the product p_n =
              // p^(n)(x)/k! into the actual form of the derivative
              {
                number k_faculty = 1;
                for (unsigned int k = 0; k <= n_derivatives; ++k)
                  {
                    values[k] *= k_faculty * lagrange_weight;
                    k_faculty *= static_cast<number>(k + 1);
                  }
              }
              break;

            // manually implement case 0 (values only), case 1 (value + first
            // derivative), and case 2 (up to second derivative) since they
            // might be called often. then, we can unroll the loop.
            case 0:
              values[0] = 1;
              for (unsigned int i = 0; i < n_supp; ++i)
                {
                  const number v = x - lagrange_support_points[i];
                  values[0] *= v;
                }
              values[0] *= lagrange_weight;
              break;

            case 1:
              values[0] = 1;
              values[1] = 0;
              for (unsigned int i = 0; i < n_supp; ++i)
                {
                  const number v = x - lagrange_support_points[i];
                  values[1]      = values[1] * v + values[0];
                  values[0] *= v;
                }
              values[0] *= lagrange_weight;
              values[1] *= lagrange_weight;
              break;

            case 2:
              values[0] = 1;
              values[1] = 0;
              values[2] = 0;
              for (unsigned int i = 0; i < n_supp; ++i)
                {
                  const number v = x - lagrange_support_points[i];
                  values[2]      = values[2] * v + values[1];
                  values[1]      = values[1] * v + values[0];
                  values[0] *= v;
                }
              values[0] *= lagrange_weight;
              values[1] *= lagrange_weight;
              values[2] *= static_cast<number>(2) * lagrange_weight;
              break;
          }
        return;
      }

    Assert(coefficients.size() > 0, ExcEmptyObject());

    // if we only need the value, then call the other function since that is
    // significantly faster (there is no need to allocate and free memory,
    // which is really expensive compared to all the other operations!)
    if (n_derivatives == 0)
      {
        values[0] = value(x);
        return;
      }

    // if there are derivatives needed, then do it properly by the full Horner
    // scheme
    const unsigned int  m = coefficients.size();
    std::vector<number> a(coefficients);
    unsigned int        j_faculty = 1;

    // loop over all requested derivatives. note that derivatives @p{j>m} are
    // necessarily zero, as they differentiate the polynomial more often than
    // the highest power is
    const unsigned int min_valuessize_m = std::min(n_derivatives + 1, m);
    for (unsigned int j = 0; j < min_valuessize_m; ++j)
      {
        for (int k = m - 2; k >= static_cast<int>(j); --k)
          a[k] += x * a[k + 1];
        values[j] = static_cast<number>(j_faculty) * a[j];

        j_faculty *= j + 1;
      }

    // fill higher derivatives by zero
    for (unsigned int j = min_valuessize_m; j <= n_derivatives; ++j)
      values[j] = 0;
  }



  template <typename number>
  void
  Polynomial<number>::transform_into_standard_form()
  {
    // should only be called when the product form is active
    Assert(in_lagrange_product_form == true, ExcInternalError());
    Assert(coefficients.size() == 0, ExcInternalError());

    // compute coefficients by expanding the product (x-x_i) term by term
    coefficients.resize(lagrange_support_points.size() + 1);
    if (lagrange_support_points.size() == 0)
      coefficients[0] = 1.;
    else
      {
        coefficients[0] = -lagrange_support_points[0];
        coefficients[1] = 1.;
        for (unsigned int i = 1; i < lagrange_support_points.size(); ++i)
          {
            coefficients[i + 1] = 1.;
            for (unsigned int j = i; j > 0; --j)
              coefficients[j] = (-lagrange_support_points[i] * coefficients[j] +
                                 coefficients[j - 1]);
            coefficients[0] *= -lagrange_support_points[i];
          }
      }
    for (unsigned int i = 0; i < lagrange_support_points.size() + 1; ++i)
      coefficients[i] *= lagrange_weight;

    // delete the product form data
    std::vector<number> new_points;
    lagrange_support_points.swap(new_points);
    in_lagrange_product_form = false;
    lagrange_weight          = 1.;
  }



  template <typename number>
  void
  Polynomial<number>::scale(std::vector<number> &coefficients,
                            const number         factor)
  {
    number f = 1.;
    for (typename std::vector<number>::iterator c = coefficients.begin();
         c != coefficients.end();
         ++c)
      {
        *c *= f;
        f *= factor;
      }
  }



  template <typename number>
  void
  Polynomial<number>::scale(const number factor)
  {
    // to scale (x-x_0)*(x-x_1)*...*(x-x_n), scale
    // support points by 1./factor and the weight
    // likewise
    if (in_lagrange_product_form == true)
      {
        number inv_fact         = number(1.) / factor;
        number accumulated_fact = 1.;
        for (unsigned int i = 0; i < lagrange_support_points.size(); ++i)
          {
            lagrange_support_points[i] *= inv_fact;
            accumulated_fact *= factor;
          }
        lagrange_weight *= accumulated_fact;
      }
    // otherwise, use the function above
    else
      scale(coefficients, factor);
  }



  template <typename number>
  void
  Polynomial<number>::multiply(std::vector<number> &coefficients,
                               const number         factor)
  {
    for (typename std::vector<number>::iterator c = coefficients.begin();
         c != coefficients.end();
         ++c)
      *c *= factor;
  }



  template <typename number>
  Polynomial<number> &
  Polynomial<number>::operator*=(const double s)
  {
    if (in_lagrange_product_form == true)
      lagrange_weight *= s;
    else
      {
        for (typename std::vector<number>::iterator c = coefficients.begin();
             c != coefficients.end();
             ++c)
          *c *= s;
      }
    return *this;
  }



  template <typename number>
  Polynomial<number> &
  Polynomial<number>::operator*=(const Polynomial<number> &p)
  {
    // if we are in Lagrange form, just append the
    // new points
    if (in_lagrange_product_form == true && p.in_lagrange_product_form == true)
      {
        lagrange_weight *= p.lagrange_weight;
        lagrange_support_points.insert(lagrange_support_points.end(),
                                       p.lagrange_support_points.begin(),
                                       p.lagrange_support_points.end());
      }

    // cannot retain product form, recompute...
    else if (in_lagrange_product_form == true)
      transform_into_standard_form();

    // need to transform p into standard form as
    // well if necessary. copy the polynomial to
    // do this
    std::unique_ptr<Polynomial<number>> q_data;
    const Polynomial<number> *          q = nullptr;
    if (p.in_lagrange_product_form == true)
      {
        q_data = std::make_unique<Polynomial<number>>(p);
        q_data->transform_into_standard_form();
        q = q_data.get();
      }
    else
      q = &p;

    // Degree of the product
    unsigned int new_degree = this->degree() + q->degree();

    std::vector<number> new_coefficients(new_degree + 1, 0.);

    for (unsigned int i = 0; i < q->coefficients.size(); ++i)
      for (unsigned int j = 0; j < this->coefficients.size(); ++j)
        new_coefficients[i + j] += this->coefficients[j] * q->coefficients[i];
    this->coefficients = new_coefficients;

    return *this;
  }



  template <typename number>
  Polynomial<number> &
  Polynomial<number>::operator+=(const Polynomial<number> &p)
  {
    // Lagrange product form cannot reasonably be
    // retained after polynomial addition. we
    // could in theory check if either this
    // polynomial or the other is a zero
    // polynomial and retain it, but we actually
    // currently (r23974) assume that the addition
    // of a zero polynomial changes the state and
    // tests equivalence.
    if (in_lagrange_product_form == true)
      transform_into_standard_form();

    // need to transform p into standard form as
    // well if necessary. copy the polynomial to
    // do this
    std::unique_ptr<Polynomial<number>> q_data;
    const Polynomial<number> *          q = nullptr;
    if (p.in_lagrange_product_form == true)
      {
        q_data = std::make_unique<Polynomial<number>>(p);
        q_data->transform_into_standard_form();
        q = q_data.get();
      }
    else
      q = &p;

    // if necessary expand the number
    // of coefficients we store
    if (q->coefficients.size() > coefficients.size())
      coefficients.resize(q->coefficients.size(), 0.);

    for (unsigned int i = 0; i < q->coefficients.size(); ++i)
      coefficients[i] += q->coefficients[i];

    return *this;
  }



  template <typename number>
  Polynomial<number> &
  Polynomial<number>::operator-=(const Polynomial<number> &p)
  {
    // Lagrange product form cannot reasonably be
    // retained after polynomial addition
    if (in_lagrange_product_form == true)
      transform_into_standard_form();

    // need to transform p into standard form as
    // well if necessary. copy the polynomial to
    // do this
    std::unique_ptr<Polynomial<number>> q_data;
    const Polynomial<number> *          q = nullptr;
    if (p.in_lagrange_product_form == true)
      {
        q_data = std::make_unique<Polynomial<number>>(p);
        q_data->transform_into_standard_form();
        q = q_data.get();
      }
    else
      q = &p;

    // if necessary expand the number
    // of coefficients we store
    if (q->coefficients.size() > coefficients.size())
      coefficients.resize(q->coefficients.size(), 0.);

    for (unsigned int i = 0; i < q->coefficients.size(); ++i)
      coefficients[i] -= q->coefficients[i];

    return *this;
  }



  template <typename number>
  bool
  Polynomial<number>::operator==(const Polynomial<number> &p) const
  {
    // need to distinguish a few cases based on
    // whether we are in product form or not. two
    // polynomials can still be the same when they
    // are on different forms, but the expansion
    // is the same
    if (in_lagrange_product_form == true && p.in_lagrange_product_form == true)
      return ((lagrange_weight == p.lagrange_weight) &&
              (lagrange_support_points == p.lagrange_support_points));
    else if (in_lagrange_product_form == true)
      {
        Polynomial<number> q = *this;
        q.transform_into_standard_form();
        return (q.coefficients == p.coefficients);
      }
    else if (p.in_lagrange_product_form == true)
      {
        Polynomial<number> q = p;
        q.transform_into_standard_form();
        return (q.coefficients == coefficients);
      }
    else
      return (p.coefficients == coefficients);
  }



  template <typename number>
  template <typename number2>
  void
  Polynomial<number>::shift(std::vector<number> &coefficients,
                            const number2        offset)
  {
    // too many coefficients cause overflow in
    // the binomial coefficient used below
    Assert(coefficients.size() < 31, ExcNotImplemented());

    // Copy coefficients to a vector of
    // accuracy given by the argument
    std::vector<number2> new_coefficients(coefficients.begin(),
                                          coefficients.end());

    // Traverse all coefficients from
    // c_1. c_0 will be modified by
    // higher degrees, only.
    for (unsigned int d = 1; d < new_coefficients.size(); ++d)
      {
        const unsigned int n = d;
        // Binomial coefficients are
        // needed for the
        // computation. The rightmost
        // value is unity.
        unsigned int binomial_coefficient = 1;

        // Powers of the offset will be
        // needed and computed
        // successively.
        number2 offset_power = offset;

        // Compute (x+offset)^d
        // and modify all values c_k
        // with k<d.
        // The coefficient in front of
        // x^d is not modified in this step.
        for (unsigned int k = 0; k < d; ++k)
          {
            // Recursion from Bronstein
            // Make sure no remainders
            // occur in integer
            // division.
            binomial_coefficient = (binomial_coefficient * (n - k)) / (k + 1);

            new_coefficients[d - k - 1] +=
              new_coefficients[d] * binomial_coefficient * offset_power;
            offset_power *= offset;
          }
        // The binomial coefficient
        // should have gone through a
        // whole row of Pascal's
        // triangle.
        Assert(binomial_coefficient == 1, ExcInternalError());
      }

    // copy new elements to old vector
    coefficients.assign(new_coefficients.begin(), new_coefficients.end());
  }



  template <typename number>
  template <typename number2>
  void
  Polynomial<number>::shift(const number2 offset)
  {
    // shift is simple for a polynomial in product
    // form, (x-x_0)*(x-x_1)*...*(x-x_n). just add
    // offset to all shifts
    if (in_lagrange_product_form == true)
      {
        for (unsigned int i = 0; i < lagrange_support_points.size(); ++i)
          lagrange_support_points[i] -= offset;
      }
    else
      // do the shift in any case
      shift(coefficients, offset);
  }



  template <typename number>
  Polynomial<number>
  Polynomial<number>::derivative() const
  {
    // no simple form possible for Lagrange
    // polynomial on product form
    if (degree() == 0)
      return Monomial<number>(0, 0.);

    std::unique_ptr<Polynomial<number>> q_data;
    const Polynomial<number> *          q = nullptr;
    if (in_lagrange_product_form == true)
      {
        q_data = std::make_unique<Polynomial<number>>(*this);
        q_data->transform_into_standard_form();
        q = q_data.get();
      }
    else
      q = this;

    std::vector<number> newcoefficients(q->coefficients.size() - 1);
    for (unsigned int i = 1; i < q->coefficients.size(); ++i)
      newcoefficients[i - 1] = number(i) * q->coefficients[i];

    return Polynomial<number>(newcoefficients);
  }



  template <typename number>
  Polynomial<number>
  Polynomial<number>::primitive() const
  {
    // no simple form possible for Lagrange
    // polynomial on product form
    std::unique_ptr<Polynomial<number>> q_data;
    const Polynomial<number> *          q = nullptr;
    if (in_lagrange_product_form == true)
      {
        q_data = std::make_unique<Polynomial<number>>(*this);
        q_data->transform_into_standard_form();
        q = q_data.get();
      }
    else
      q = this;

    std::vector<number> newcoefficients(q->coefficients.size() + 1);
    newcoefficients[0] = 0.;
    for (unsigned int i = 0; i < q->coefficients.size(); ++i)
      newcoefficients[i + 1] = q->coefficients[i] / number(i + 1.);

    return Polynomial<number>(newcoefficients);
  }



  template <typename number>
  void
  Polynomial<number>::print(std::ostream &out) const
  {
    if (in_lagrange_product_form == true)
      {
        out << lagrange_weight;
        for (unsigned int i = 0; i < lagrange_support_points.size(); ++i)
          out << " (x-" << lagrange_support_points[i] << ")";
        out << std::endl;
      }
    else
      for (int i = degree(); i >= 0; --i)
        {
          out << coefficients[i] << " x^" << i << std::endl;
        }
  }


  template <typename number>
  std::size_t
  Polynomial<number>::memory_consumption() const
  {
    return (MemoryConsumption::memory_consumption(coefficients) +
            MemoryConsumption::memory_consumption(in_lagrange_product_form) +
            MemoryConsumption::memory_consumption(lagrange_support_points) +
            MemoryConsumption::memory_consumption(lagrange_weight));
  }



  // ------------------ class Monomial -------------------------- //

  template <typename number>
  std::vector<number>
  Monomial<number>::make_vector(unsigned int n, double coefficient)
  {
    std::vector<number> result(n + 1, 0.);
    result[n] = coefficient;
    return result;
  }



  template <typename number>
  Monomial<number>::Monomial(unsigned int n, double coefficient)
    : Polynomial<number>(make_vector(n, coefficient))
  {}



  template <typename number>
  std::vector<Polynomial<number>>
  Monomial<number>::generate_complete_basis(const unsigned int degree)
  {
    std::vector<Polynomial<number>> v;
    v.reserve(degree + 1);
    for (unsigned int i = 0; i <= degree; ++i)
      v.push_back(Monomial<number>(i));
    return v;
  }



  // ------------------ class LagrangeEquidistant --------------- //

  namespace internal
  {
    namespace LagrangeEquidistantImplementation
    {
      std::vector<Point<1>>
      generate_equidistant_unit_points(const unsigned int n)
      {
        std::vector<Point<1>> points(n + 1);
        const double          one_over_n = 1. / n;
        for (unsigned int k = 0; k <= n; ++k)
          points[k](0) = static_cast<double>(k) * one_over_n;
        return points;
      }
    } // namespace LagrangeEquidistantImplementation
  }   // namespace internal



  LagrangeEquidistant::LagrangeEquidistant(const unsigned int n,
                                           const unsigned int support_point)
    : Polynomial<double>(internal::LagrangeEquidistantImplementation::
                           generate_equidistant_unit_points(n),
                         support_point)
  {
    Assert(coefficients.size() == 0, ExcInternalError());

    // For polynomial order up to 3, we have precomputed weights. Use these
    // weights instead of the product form
    if (n <= 3)
      {
        this->in_lagrange_product_form = false;
        this->lagrange_weight          = 1.;
        std::vector<double> new_support_points;
        this->lagrange_support_points.swap(new_support_points);
        this->coefficients.resize(n + 1);
        compute_coefficients(n, support_point, this->coefficients);
      }
  }



  void
  LagrangeEquidistant::compute_coefficients(const unsigned int   n,
                                            const unsigned int   support_point,
                                            std::vector<double> &a)
  {
    AssertIndexRange(support_point, n + 1);

    unsigned int n_functions = n + 1;
    AssertIndexRange(support_point, n_functions);
    double const *x = nullptr;

    switch (n)
      {
        case 1:
          {
            static const double x1[4] = {1.0, -1.0, 0.0, 1.0};
            x                         = &x1[0];
            break;
          }
        case 2:
          {
            static const double x2[9] = {
              1.0, -3.0, 2.0, 0.0, 4.0, -4.0, 0.0, -1.0, 2.0};
            x = &x2[0];
            break;
          }
        case 3:
          {
            static const double x3[16] = {1.0,
                                          -11.0 / 2.0,
                                          9.0,
                                          -9.0 / 2.0,
                                          0.0,
                                          9.0,
                                          -45.0 / 2.0,
                                          27.0 / 2.0,
                                          0.0,
                                          -9.0 / 2.0,
                                          18.0,
                                          -27.0 / 2.0,
                                          0.0,
                                          1.0,
                                          -9.0 / 2.0,
                                          9.0 / 2.0};
            x                          = &x3[0];
            break;
          }
        default:
          Assert(false, ExcInternalError())
      }

    Assert(x != nullptr, ExcInternalError());
    for (unsigned int i = 0; i < n_functions; ++i)
      a[i] = x[support_point * n_functions + i];
  }



  std::vector<Polynomial<double>>
  LagrangeEquidistant::generate_complete_basis(const unsigned int degree)
  {
    if (degree == 0)
      // create constant polynomial
      return std::vector<Polynomial<double>>(
        1, Polynomial<double>(std::vector<double>(1, 1.)));
    else
      {
        // create array of Lagrange
        // polynomials
        std::vector<Polynomial<double>> v;
        for (unsigned int i = 0; i <= degree; ++i)
          v.push_back(LagrangeEquidistant(degree, i));
        return v;
      }
  }



  //----------------------------------------------------------------------//


  std::vector<Polynomial<double>>
  generate_complete_Lagrange_basis(const std::vector<Point<1>> &points)
  {
    std::vector<Polynomial<double>> p;
    p.reserve(points.size());

    for (unsigned int i = 0; i < points.size(); ++i)
      p.emplace_back(points, i);
    return p;
  }



  // ------------------ class Legendre --------------- //



  Legendre::Legendre(const unsigned int k)
    : Polynomial<double>(0)
  {
    this->coefficients.clear();
    this->in_lagrange_product_form = true;
    this->lagrange_support_points.resize(k);

    // the roots of a Legendre polynomial are exactly the points in the
    // Gauss-Legendre quadrature formula
    if (k > 0)
      {
        QGauss<1> gauss(k);
        for (unsigned int i = 0; i < k; ++i)
          this->lagrange_support_points[i] = gauss.get_points()[i][0];
      }

    // compute the abscissa in zero of the product of monomials. The exact
    // value should be sqrt(2*k+1), so set the weight to that value.
    double prod = 1.;
    for (unsigned int i = 0; i < k; ++i)
      prod *= this->lagrange_support_points[i];
    this->lagrange_weight = std::sqrt(double(2 * k + 1)) / prod;
  }



  std::vector<Polynomial<double>>
  Legendre::generate_complete_basis(const unsigned int degree)
  {
    std::vector<Polynomial<double>> v;
    v.reserve(degree + 1);
    for (unsigned int i = 0; i <= degree; ++i)
      v.push_back(Legendre(i));
    return v;
  }



  // ------------------ class Lobatto -------------------- //


  Lobatto::Lobatto(const unsigned int p)
    : Polynomial<double>(compute_coefficients(p))
  {}

  std::vector<double>
  Lobatto::compute_coefficients(const unsigned int p)
  {
    switch (p)
      {
        case 0:
          {
            std::vector<double> coefficients(2);

            coefficients[0] = 1.0;
            coefficients[1] = -1.0;
            return coefficients;
          }

        case 1:
          {
            std::vector<double> coefficients(2);

            coefficients[0] = 0.0;
            coefficients[1] = 1.0;
            return coefficients;
          }

        case 2:
          {
            std::vector<double> coefficients(3);

            coefficients[0] = 0.0;
            coefficients[1] = -1.0 * std::sqrt(3.);
            coefficients[2] = std::sqrt(3.);
            return coefficients;
          }

        default:
          {
            std::vector<double> coefficients(p + 1);
            std::vector<double> legendre_coefficients_tmp1(p);
            std::vector<double> legendre_coefficients_tmp2(p - 1);

            coefficients[0]               = -1.0 * std::sqrt(3.);
            coefficients[1]               = 2.0 * std::sqrt(3.);
            legendre_coefficients_tmp1[0] = 1.0;

            for (unsigned int i = 2; i < p; ++i)
              {
                for (unsigned int j = 0; j < i - 1; ++j)
                  legendre_coefficients_tmp2[j] = legendre_coefficients_tmp1[j];

                for (unsigned int j = 0; j < i; ++j)
                  legendre_coefficients_tmp1[j] = coefficients[j];

                coefficients[0] =
                  std::sqrt(2 * i + 1.) *
                  ((1.0 - 2 * i) * legendre_coefficients_tmp1[0] /
                     std::sqrt(2 * i - 1.) +
                   (1.0 - i) * legendre_coefficients_tmp2[0] /
                     std::sqrt(2 * i - 3.)) /
                  i;

                for (unsigned int j = 1; j < i - 1; ++j)
                  coefficients[j] =
                    std::sqrt(2 * i + 1.) *
                    (std::sqrt(2 * i - 1.) *
                       (2.0 * legendre_coefficients_tmp1[j - 1] -
                        legendre_coefficients_tmp1[j]) +
                     (1.0 - i) * legendre_coefficients_tmp2[j] /
                       std::sqrt(2 * i - 3.)) /
                    i;

                coefficients[i - 1] = std::sqrt(4 * i * i - 1.) *
                                      (2.0 * legendre_coefficients_tmp1[i - 2] -
                                       legendre_coefficients_tmp1[i - 1]) /
                                      i;
                coefficients[i] = 2.0 * std::sqrt(4 * i * i - 1.) *
                                  legendre_coefficients_tmp1[i - 1] / i;
              }

            for (int i = p; i > 0; --i)
              coefficients[i] = coefficients[i - 1] / i;

            coefficients[0] = 0.0;
            return coefficients;
          }
      }
  }

  std::vector<Polynomial<double>>
  Lobatto::generate_complete_basis(const unsigned int p)
  {
    std::vector<Polynomial<double>> basis(p + 1);

    for (unsigned int i = 0; i <= p; ++i)
      basis[i] = Lobatto(i);

    return basis;
  }



  // ------------------ class Hierarchical --------------- //


  // Reserve space for polynomials up to degree 19. Should be sufficient
  // for the start.
  std::vector<std::unique_ptr<const std::vector<double>>>
    Hierarchical::recursive_coefficients(20);



  Hierarchical::Hierarchical(const unsigned int k)
    : Polynomial<double>(get_coefficients(k))
  {}



  void
  Hierarchical::compute_coefficients(const unsigned int k_)
  {
    unsigned int k = k_;

    // first make sure that no other
    // thread intercepts the operation
    // of this function
    // for this, acquire the lock
    // until we quit this function
    std::lock_guard<std::mutex> lock(coefficients_lock);

    // The first 2 coefficients
    // are hard-coded
    if (k == 0)
      k = 1;
    // check: does the information
    // already exist?
    if ((recursive_coefficients.size() < k + 1) ||
        (recursive_coefficients[k].get() == nullptr))
      // no, then generate the
      // respective coefficients
      {
        // make sure that there is enough
        // space in the array for the
        // coefficients, so we have to resize
        // it to size k+1

        // but it's more complicated than
        // that: we call this function
        // recursively, so if we simply
        // resize it to k+1 here, then
        // compute the coefficients for
        // degree k-1 by calling this
        // function recursively, then it will
        // reset the size to k -- not enough
        // for what we want to do below. the
        // solution therefore is to only
        // resize the size if we are going to
        // *increase* it
        if (recursive_coefficients.size() < k + 1)
          recursive_coefficients.resize(k + 1);

        if (k <= 1)
          {
            // create coefficients
            // vectors for k=0 and k=1
            //
            // allocate the respective
            // amount of memory and
            // later assign it to the
            // coefficients array to
            // make it const
            std::vector<double> c0(2);
            c0[0] = 1.;
            c0[1] = -1.;

            std::vector<double> c1(2);
            c1[0] = 0.;
            c1[1] = 1.;

            // now make these arrays
            // const
            recursive_coefficients[0] =
              std::make_unique<const std::vector<double>>(std::move(c0));
            recursive_coefficients[1] =
              std::make_unique<const std::vector<double>>(std::move(c1));
          }
        else if (k == 2)
          {
            coefficients_lock.unlock();
            compute_coefficients(1);
            coefficients_lock.lock();

            std::vector<double> c2(3);

            const double a = 1.; // 1./8.;

            c2[0] = 0. * a;
            c2[1] = -4. * a;
            c2[2] = 4. * a;

            recursive_coefficients[2] =
              std::make_unique<const std::vector<double>>(std::move(c2));
          }
        else
          {
            // for larger numbers,
            // compute the coefficients
            // recursively. to do so,
            // we have to release the
            // lock temporarily to
            // allow the called
            // function to acquire it
            // itself
            coefficients_lock.unlock();
            compute_coefficients(k - 1);
            coefficients_lock.lock();

            std::vector<double> ck(k + 1);

            const double a = 1.; // 1./(2.*k);

            ck[0] = -a * (*recursive_coefficients[k - 1])[0];

            for (unsigned int i = 1; i <= k - 1; ++i)
              ck[i] = a * (2. * (*recursive_coefficients[k - 1])[i - 1] -
                           (*recursive_coefficients[k - 1])[i]);

            ck[k] = a * 2. * (*recursive_coefficients[k - 1])[k - 1];
            // for even degrees, we need
            // to add a multiple of
            // basis fcn phi_2
            if ((k % 2) == 0)
              {
                double b = 1.; // 8.;
                // for (unsigned int i=1; i<=k; i++)
                //  b /= 2.*i;

                ck[1] += b * (*recursive_coefficients[2])[1];
                ck[2] += b * (*recursive_coefficients[2])[2];
              }
            // finally assign the newly
            // created vector to the
            // const pointer in the
            // coefficients array
            recursive_coefficients[k] =
              std::make_unique<const std::vector<double>>(std::move(ck));
          }
      }
  }



  const std::vector<double> &
  Hierarchical::get_coefficients(const unsigned int k)
  {
    // first make sure the coefficients
    // get computed if so necessary
    compute_coefficients(k);

    // then get a pointer to the array
    // of coefficients. do that in a MT
    // safe way
    std::lock_guard<std::mutex> lock(coefficients_lock);
    return *recursive_coefficients[k];
  }



  std::vector<Polynomial<double>>
  Hierarchical::generate_complete_basis(const unsigned int degree)
  {
    if (degree == 0)
      // create constant
      // polynomial. note that we
      // can't use the other branch
      // of the if-statement, since
      // calling the constructor of
      // this class with argument
      // zero does _not_ create the
      // constant polynomial, but
      // rather 1-x
      return std::vector<Polynomial<double>>(
        1, Polynomial<double>(std::vector<double>(1, 1.)));
    else
      {
        std::vector<Polynomial<double>> v;
        v.reserve(degree + 1);
        for (unsigned int i = 0; i <= degree; ++i)
          v.push_back(Hierarchical(i));
        return v;
      }
  }



  // ------------------ HermiteInterpolation --------------- //

  HermiteInterpolation::HermiteInterpolation(const unsigned int p)
    : Polynomial<double>(0)
  {
    this->coefficients.clear();
    this->in_lagrange_product_form = true;

    this->lagrange_support_points.resize(3);
    if (p == 0)
      {
        this->lagrange_support_points[0] = -0.5;
        this->lagrange_support_points[1] = 1.;
        this->lagrange_support_points[2] = 1.;
        this->lagrange_weight            = 2.;
      }
    else if (p == 1)
      {
        this->lagrange_support_points[0] = 0.;
        this->lagrange_support_points[1] = 0.;
        this->lagrange_support_points[2] = 1.5;
        this->lagrange_weight            = -2.;
      }
    else if (p == 2)
      {
        this->lagrange_support_points[0] = 0.;
        this->lagrange_support_points[1] = 1.;
        this->lagrange_support_points[2] = 1.;
      }
    else if (p == 3)
      {
        this->lagrange_support_points[0] = 0.;
        this->lagrange_support_points[1] = 0.;
        this->lagrange_support_points[2] = 1.;
      }
    else
      {
        this->lagrange_support_points.resize(4);
        this->lagrange_support_points[0] = 0.;
        this->lagrange_support_points[1] = 0.;
        this->lagrange_support_points[2] = 1.;
        this->lagrange_support_points[3] = 1.;
        this->lagrange_weight            = 16.;

        if (p > 4)
          {
            Legendre legendre(p - 4);
            (*this) *= legendre;
          }
      }
  }


  std::vector<Polynomial<double>>
  HermiteInterpolation::generate_complete_basis(const unsigned int n)
  {
    Assert(n >= 3,
           ExcNotImplemented("Hermite interpolation makes no sense for "
                             "degrees less than three"));
    std::vector<Polynomial<double>> basis(n + 1);

    for (unsigned int i = 0; i <= n; ++i)
      basis[i] = HermiteInterpolation(i);

    return basis;
  }


  // ------------------ HermiteLikeInterpolation --------------- //
  namespace
  {
    // Finds the zero position x_star such that the mass matrix entry (0,1)
    // with the Hermite polynomials evaluates to zero. The function has
    // originally been derived by a secant method for the integral entry
    // l_0(x) * l_1(x) but we only need to do one iteration because the zero
    // x_star is linear in the integral value.
    double
    find_support_point_x_star(const std::vector<double> &jacobi_roots)
    {
      // Initial guess for the support point position values: The zero turns
      // out to be between zero and the first root of the Jacobi polynomial,
      // but the algorithm is agnostic about that, so simply choose two points
      // that are sufficiently far apart.
      double             guess_left  = 0;
      double             guess_right = 0.5;
      const unsigned int degree      = jacobi_roots.size() + 3;

      // Compute two integrals of the product of l_0(x) * l_1(x)
      // l_0(x) =
      // (x-y)*(x-jacobi_roots(0))*...*(x-jacobi_roos(degree-4))*(x-1)*(x-1)
      // l_1(x) =
      // (x-0)*(x-jacobi_roots(0))*...*(x-jacobi_roots(degree-4))*(x-1)*(x-1)
      // where y is either guess_left or guess_right for the two integrals.
      // Note that the polynomials are not yet normalized here, which is not
      // necessary because we are only looking for the x_star where the matrix
      // entry is zero, for which the constants do not matter.
      QGauss<1> gauss(degree + 1);
      double    integral_left = 0, integral_right = 0;
      for (unsigned int q = 0; q < gauss.size(); ++q)
        {
          const double x               = gauss.point(q)[0];
          double       poly_val_common = x;
          for (unsigned int j = 0; j < degree - 3; ++j)
            poly_val_common *= Utilities::fixed_power<2>(x - jacobi_roots[j]);
          poly_val_common *= Utilities::fixed_power<4>(x - 1.);
          integral_left +=
            gauss.weight(q) * (poly_val_common * (x - guess_left));
          integral_right +=
            gauss.weight(q) * (poly_val_common * (x - guess_right));
        }

      // compute guess by secant method. Due to linearity in the root x_star,
      // this is the correct position after this single step
      return guess_right - (guess_right - guess_left) /
                             (integral_right - integral_left) * integral_right;
    }
  } // namespace



  HermiteLikeInterpolation::HermiteLikeInterpolation(const unsigned int degree,
                                                     const unsigned int index)
    : Polynomial<double>(0)
  {
    AssertIndexRange(index, degree + 1);

    this->coefficients.clear();
    this->in_lagrange_product_form = true;

    this->lagrange_support_points.resize(degree);

    if (degree == 0)
      this->lagrange_weight = 1.;
    else if (degree == 1)
      {
        if (index == 0)
          {
            this->lagrange_support_points[0] = 1.;
            this->lagrange_weight            = -1.;
          }
        else
          {
            this->lagrange_support_points[0] = 0.;
            this->lagrange_weight            = 1.;
          }
      }
    else if (degree == 2)
      {
        if (index == 0)
          {
            this->lagrange_support_points[0] = 1.;
            this->lagrange_support_points[1] = 1.;
            this->lagrange_weight            = 1.;
          }
        else if (index == 1)
          {
            this->lagrange_support_points[0] = 0;
            this->lagrange_support_points[1] = 1;
            this->lagrange_weight            = -2.;
          }
        else
          {
            this->lagrange_support_points[0] = 0.;
            this->lagrange_support_points[1] = 0.;
            this->lagrange_weight            = 1.;
          }
      }
    else if (degree == 3)
      {
        // 4 Polynomials with degree 3
        // entries (1,0) and (3,2) of the mass matrix will be equal to 0
        //
        //     | x  0  x  x |
        //     | 0  x  x  x |
        // M = | x  x  x  0 |
        //     | x  x  0  x |
        //
        if (index == 0)
          {
            this->lagrange_support_points[0] = 2. / 7.;
            this->lagrange_support_points[1] = 1.;
            this->lagrange_support_points[2] = 1.;
            this->lagrange_weight            = -3.5;
          }
        else if (index == 1)
          {
            this->lagrange_support_points[0] = 0.;
            this->lagrange_support_points[1] = 1.;
            this->lagrange_support_points[2] = 1.;

            // this magic value 5.5 is obtained when evaluating the general
            // formula below for the degree=3 case
            this->lagrange_weight = 5.5;
          }
        else if (index == 2)
          {
            this->lagrange_support_points[0] = 0.;
            this->lagrange_support_points[1] = 0.;
            this->lagrange_support_points[2] = 1.;
            this->lagrange_weight            = -5.5;
          }
        else if (index == 3)
          {
            this->lagrange_support_points[0] = 0.;
            this->lagrange_support_points[1] = 0.;
            this->lagrange_support_points[2] = 5. / 7.;
            this->lagrange_weight            = 3.5;
          }
      }
    else
      {
        // Higher order Polynomials degree>=4: the entries (1,0) and
        // (degree,degree-1) of the mass matrix will be equal to 0
        //
        //     | x  0  x  x         x  x  x |
        //     | 0  x  x  x  . . .  x  x  x |
        //     | x  x  x  0         0  x  x |
        //     | x  x  0  x         0  x  x |
        //     |     .       .         .    |
        // M = |     .         .       .    |
        //     |     .           .     .    |
        //     | x  x  0  0         x  x  x |
        //     | x  x  x  x  . . .  x  x  0 |
        //     | x  x  x  x         x  0  x |
        //
        // We find the inner points as the zeros of the Jacobi polynomials
        // with alpha = beta = 4 which is the polynomial with the kernel
        // (1-x)^4 (1+x)^4. Since polynomials (1-x)^2 (1+x)^2 are contained
        // in every interior polynomial (bubble function), their product
        // leads us to the orthogonality condition of the Jacobi(4,4)
        // polynomials.

        std::vector<double> jacobi_roots =
          jacobi_polynomial_roots<double>(degree - 3, 4, 4);
        AssertDimension(jacobi_roots.size(), degree - 3);

        // iteration from variable support point N with secant method
        // initial values

        this->lagrange_support_points.resize(degree);
        if (index == 0)
          {
            const double auxiliary_zero =
              find_support_point_x_star(jacobi_roots);
            this->lagrange_support_points[0] = auxiliary_zero;
            for (unsigned int m = 0; m < degree - 3; m++)
              this->lagrange_support_points[m + 1] = jacobi_roots[m];
            this->lagrange_support_points[degree - 2] = 1.;
            this->lagrange_support_points[degree - 1] = 1.;

            // ensure that the polynomial evaluates to one at x=0
            this->lagrange_weight = 1. / this->value(0);
          }
        else if (index == 1)
          {
            this->lagrange_support_points[0] = 0.;
            for (unsigned int m = 0; m < degree - 3; m++)
              this->lagrange_support_points[m + 1] = jacobi_roots[m];
            this->lagrange_support_points[degree - 2] = 1.;
            this->lagrange_support_points[degree - 1] = 1.;

            // Select the weight to make the derivative of the sum of P_0 and
            // P_1 in zero to be 0. The derivative in x=0 is simply given by
            // p~(0)/auxiliary_zero+p~'(0) + a*p~(0), where p~(x) is the
            // Lagrange polynomial in all points except the first one which is
            // the same for P_0 and P_1, and a is the weight we seek here. If
            // we solve this for a, we obtain the desired property. Since the
            // basis is nodal for all interior points, this property ensures
            // that the sum of all polynomials with weight 1 is one.
            std::vector<Point<1>> points(degree);
            double                ratio = 1.;
            for (unsigned int i = 0; i < degree; ++i)
              {
                points[i][0] = this->lagrange_support_points[i];
                if (i > 0)
                  ratio *= -this->lagrange_support_points[i];
              }
            Polynomial<double>  helper(points, 0);
            std::vector<double> value_and_grad(2);
            helper.value(0., value_and_grad);
            Assert(std::abs(value_and_grad[0]) > 1e-10,
                   ExcInternalError("There should not be a zero at x=0."));

            const double auxiliary_zero =
              find_support_point_x_star(jacobi_roots);
            this->lagrange_weight =
              (1. / auxiliary_zero - value_and_grad[1] / value_and_grad[0]) /
              ratio;
          }
        else if (index >= 2 && index < degree - 1)
          {
            this->lagrange_support_points[0] = 0.;
            this->lagrange_support_points[1] = 0.;
            for (unsigned int m = 0, c = 2; m < degree - 3; m++)
              if (m + 2 != index)
                this->lagrange_support_points[c++] = jacobi_roots[m];
            this->lagrange_support_points[degree - 2] = 1.;
            this->lagrange_support_points[degree - 1] = 1.;

            // ensure that the polynomial evaluates to one at the respective
            // nodal point
            this->lagrange_weight = 1. / this->value(jacobi_roots[index - 2]);
          }
        else if (index == degree - 1)
          {
            this->lagrange_support_points[0] = 0.;
            this->lagrange_support_points[1] = 0.;
            for (unsigned int m = 0; m < degree - 3; m++)
              this->lagrange_support_points[m + 2] = jacobi_roots[m];
            this->lagrange_support_points[degree - 1] = 1.;

            std::vector<Point<1>> points(degree);
            double                ratio = 1.;
            for (unsigned int i = 0; i < degree; ++i)
              {
                points[i][0] = this->lagrange_support_points[i];
                if (i < degree - 1)
                  ratio *= 1. - this->lagrange_support_points[i];
              }
            Polynomial<double>  helper(points, degree - 1);
            std::vector<double> value_and_grad(2);
            helper.value(1., value_and_grad);
            Assert(std::abs(value_and_grad[0]) > 1e-10,
                   ExcInternalError("There should not be a zero at x=1."));

            const double auxiliary_zero =
              find_support_point_x_star(jacobi_roots);
            this->lagrange_weight =
              (-1. / auxiliary_zero - value_and_grad[1] / value_and_grad[0]) /
              ratio;
          }
        else if (index == degree)
          {
            const double auxiliary_zero =
              find_support_point_x_star(jacobi_roots);
            this->lagrange_support_points[0] = 0.;
            this->lagrange_support_points[1] = 0.;
            for (unsigned int m = 0; m < degree - 3; m++)
              this->lagrange_support_points[m + 2] = jacobi_roots[m];
            this->lagrange_support_points[degree - 1] = 1. - auxiliary_zero;

            // ensure that the polynomial evaluates to one at x=1
            this->lagrange_weight = 1. / this->value(1.);
          }
      }
  }



  std::vector<Polynomial<double>>
  HermiteLikeInterpolation::generate_complete_basis(const unsigned int degree)
  {
    std::vector<Polynomial<double>> basis(degree + 1);

    for (unsigned int i = 0; i <= degree; ++i)
      basis[i] = HermiteLikeInterpolation(degree, i);

    return basis;
  }

} // namespace Polynomials

// ------------------ explicit instantiations --------------- //

#ifndef DOXYGEN
namespace Polynomials
{
  template class Polynomial<float>;
  template class Polynomial<double>;
  template class Polynomial<long double>;

  template void
  Polynomial<float>::shift(const float offset);
  template void
  Polynomial<float>::shift(const double offset);
  template void
  Polynomial<double>::shift(const double offset);
  template void
  Polynomial<long double>::shift(const long double offset);
  template void
  Polynomial<float>::shift(const long double offset);
  template void
  Polynomial<double>::shift(const long double offset);

  template class Monomial<float>;
  template class Monomial<double>;
  template class Monomial<long double>;
} // namespace Polynomials
#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE
