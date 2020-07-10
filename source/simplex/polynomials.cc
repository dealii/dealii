// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


#include <deal.II/simplex/polynomials.h>

DEAL_II_NAMESPACE_OPEN

namespace Simplex
{
  namespace
  {
    unsigned int
    compute_n_polynomials(const unsigned int dim, const unsigned int degree)
    {
      if (dim == 1)
        {
          if (degree == 1)
            return 2;
          if (degree == 2)
            return 3;
        }
      else if (dim == 2)
        {
          if (degree == 1)
            return 3;
          if (degree == 2)
            return 6;
        }
      else if (dim == 3)
        {
          if (degree == 1)
            return 4;
          if (degree == 2)
            return 10;
        }

      Assert(false, ExcNotImplemented());

      return 0;
    }
  } // namespace



  template <int dim>
  ScalarPolynomial<dim>::ScalarPolynomial(const unsigned int degree)
    : ScalarPolynomialsBase<dim>(degree, compute_n_polynomials(dim, degree))
  {}



  template <int dim>
  double
  ScalarPolynomial<dim>::compute_value(const unsigned int i,
                                       const Point<dim> & p) const
  {
    if (dim == 1)
      {
        if (this->degree() == 1)
          {
            if (i == 0)
              return 1.0 - p[0];
            else if (i == 1)
              return p[0];
          }
        else if (this->degree() == 2)
          {
            if (i == 0)
              return 2.0 * p[0] * p[0] - 3.0 * p[0] + 1;
            else if (i == 1)
              return 2.0 * p[0] * p[0] - p[0];
            else if (i == 2)
              return -4.0 * p[0] * p[0] + 4.0 * p[0];
          }
      }
    else if (dim == 2)
      {
        if (this->degree() == 1)
          {
            if (i == 0)
              return 1.0 - p[0] - p[1];
            else if (i == 1)
              return p[0];
            else if (i == 2)
              return p[1];
          }
        else if (this->degree() == 2)
          {
            const double t1 = 1.0 - p[0] - p[1];
            const double t2 = p[0];
            const double t3 = p[1];

            if (i == 0)
              return t1 * (2.0 * t1 - 1.0);
            else if (i == 1)
              return t2 * (2.0 * t2 - 1.0);
            else if (i == 2)
              return t3 * (2.0 * t3 - 1.0);
            else if (i == 3)
              return 4.0 * t2 * t1;
            else if (i == 4)
              return 4.0 * t2 * t3;
            else if (i == 5)
              return 4.0 * t3 * t1;
          }
      }
    else if (dim == 3)
      {
        if (this->degree() == 1)
          {
            if (i == 0)
              return 1.0 - p[0] - p[1] - p[2];
            else if (i == 1)
              return p[0];
            else if (i == 2)
              return p[1];
            else if (i == 3)
              return p[2];
          }
        else if (this->degree() == 2)
          {
            const double r = p[0];
            const double s = p[1];
            const double t = p[2];
            const double u = 1.0 - p[0] - p[1] - p[2];
            if (i == 0)
              return u * (2.0 * u - 1.0);
            else if (i == 1)
              return r * (2.0 * r - 1.0);
            else if (i == 2)
              return s * (2.0 * s - 1.0);
            else if (i == 3)
              return t * (2.0 * t - 1.0);
            else if (i == 4)
              return 4.0 * r * u;
            else if (i == 5)
              return 4.0 * r * s;
            else if (i == 6)
              return 4.0 * s * u;
            else if (i == 7)
              return 4.0 * t * u;
            else if (i == 8)
              return 4.0 * r * t;
            else if (i == 9)
              return 4.0 * s * t;
          }
      }

    Assert(false, ExcNotImplemented());

    return 0;
  }



  template <int dim>
  Tensor<1, dim>
  ScalarPolynomial<dim>::compute_grad(const unsigned int i,
                                      const Point<dim> & p) const
  {
    Tensor<1, dim> grad;

    if (dim == 1)
      {
        if (this->degree() == 1)
          {
            if (i == 0)
              grad[0] = -1.0;
            else if (i == 1)
              grad[0] = 1.0;
          }
        else if (this->degree() == 2)
          {
            if (i == 0)
              grad[0] = 4.0 * p[0] - 3.0;
            else if (i == 1)
              grad[0] = 4.0 * p[0] - 1.0;
            else if (i == 2)
              grad[0] = -8.0 * p[0] + 4.0;
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
    else if (dim == 2)
      {
        if (this->degree() == 1)
          {
            if (i == 0)
              {
                grad[0] = -1.0;
                grad[1] = -1.0;
              }
            else if (i == 1)
              {
                grad[0] = +1.0;
                grad[1] = +0.0;
              }
            else if (i == 2)
              {
                grad[0] = +0.0;
                grad[1] = +1.0;
              }
            else
              {
                Assert(false, ExcNotImplemented());
              }
          }
        else if (this->degree() == 2)
          {
            if (i == 0)
              {
                grad[0] = -3.0 + 4.0 * (p[0] + p[1]);
                grad[1] = -3.0 + 4.0 * (p[0] + p[1]);
              }
            else if (i == 1)
              {
                grad[0] = 4.0 * p[0] - 1.0;
                grad[1] = 0.0;
              }
            else if (i == 2)
              {
                grad[0] = 0.0;
                grad[1] = 4.0 * p[1] - 1.0;
              }
            else if (i == 3)
              {
                grad[0] = 4.0 * (1.0 - 2.0 * p[0] - p[1]);
                grad[1] = -4.0 * p[0];
              }
            else if (i == 4)
              {
                grad[0] = 4.0 * p[1];
                grad[1] = 4.0 * p[0];
              }
            else if (i == 5)
              {
                grad[0] = -4.0 * p[1];
                grad[1] = 4.0 * (1.0 - p[0] - 2.0 * p[1]);
              }
            else
              {
                Assert(false, ExcNotImplemented());
              }
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
    else if (dim == 3)
      {
        if (this->degree() == 1)
          {
            if (i == 0)
              {
                grad[0] = -1.0;
                grad[1] = -1.0;
                grad[2] = -1.0;
              }
            else if (i == 1)
              {
                grad[0] = +1.0;
                grad[1] = +0.0;
                grad[2] = +0.0;
              }
            else if (i == 2)
              {
                grad[0] = +0.0;
                grad[1] = +1.0;
                grad[2] = +0.0;
              }
            else if (i == 3)
              {
                grad[0] = +0.0;
                grad[1] = +0.0;
                grad[2] = +1.0;
              }
          }
        else if (this->degree() == 2)
          {
            const double r = p[0];
            const double s = p[1];
            const double t = p[2];
            const double u = 1.0 - p[0] - p[1] - p[2];

            if (i == 0)
              {
                grad[0] = -4.0 * u + 1.;
                grad[1] = grad[0];
                grad[2] = grad[0];
              }
            else if (i == 1)
              {
                grad[0] = +4.0 * r - 1.;
                grad[1] = +0.0;
                grad[2] = +0.0;
              }
            else if (i == 2)
              {
                grad[0] = +0.0;
                grad[1] = +4.0 * s - 1.;
                grad[2] = +0.0;
              }
            else if (i == 3)
              {
                grad[0] = +0.0;
                grad[1] = +0.0;
                grad[2] = +4.0 * t - 1.;
              }
            else if (i == 4)
              {
                grad[0] = +4.0 * (u - r);
                grad[1] = -4.0 * r;
                grad[2] = -4.0 * r;
              }
            else if (i == 5)
              {
                grad[0] = +4.0 * s;
                grad[1] = +4.0 * r;
                grad[2] = +0.0;
              }
            else if (i == 6)
              {
                grad[0] = -4.0 * s;
                grad[1] = +4.0 * (u - s);
                grad[2] = -4.0 * s;
              }
            else if (i == 7)
              {
                grad[0] = -4.0 * t;
                grad[1] = -4.0 * t;
                grad[2] = +4.0 * (u - t);
              }
            else if (i == 8)
              {
                grad[0] = +4.0 * t;
                grad[1] = +0.0;
                grad[2] = +4.0 * r;
              }
            else if (i == 9)
              {
                grad[0] = +0.0;
                grad[1] = +4.0 * t;
                grad[2] = +4.0 * s;
              }
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
    else
      {
        Assert(false, ExcNotImplemented());
      }

    return grad;
  }



  template <int dim>
  Tensor<2, dim>
  ScalarPolynomial<dim>::compute_grad_grad(const unsigned int i,
                                           const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());
    return Tensor<2, dim>();
  }



  template <int dim>
  void
  ScalarPolynomial<dim>::evaluate(
    const Point<dim> &           unit_point,
    std::vector<double> &        values,
    std::vector<Tensor<1, dim>> &grads,
    std::vector<Tensor<2, dim>> &grad_grads,
    std::vector<Tensor<3, dim>> &third_derivatives,
    std::vector<Tensor<4, dim>> &fourth_derivatives) const
  {
    (void)grads;
    (void)grad_grads;
    (void)third_derivatives;
    (void)fourth_derivatives;

    if (values.size() == this->n())
      for (unsigned int i = 0; i < this->n(); i++)
        values[i] = compute_value(i, unit_point);

    if (grads.size() == this->n())
      for (unsigned int i = 0; i < this->n(); i++)
        grads[i] = compute_grad(i, unit_point);
  }



  template <int dim>
  Tensor<1, dim>
  ScalarPolynomial<dim>::compute_1st_derivative(const unsigned int i,
                                                const Point<dim> & p) const
  {
    return compute_grad(i, p);
  }



  template <int dim>
  Tensor<2, dim>
  ScalarPolynomial<dim>::compute_2nd_derivative(const unsigned int i,
                                                const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());

    return {};
  }



  template <int dim>
  Tensor<3, dim>
  ScalarPolynomial<dim>::compute_3rd_derivative(const unsigned int i,
                                                const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());

    return {};
  }



  template <int dim>
  Tensor<4, dim>
  ScalarPolynomial<dim>::compute_4th_derivative(const unsigned int i,
                                                const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());

    return {};
  }



  template <int dim>
  std::string
  ScalarPolynomial<dim>::name() const
  {
    return "Simplex";
  }



  template <int dim>
  std::unique_ptr<ScalarPolynomialsBase<dim>>
  ScalarPolynomial<dim>::clone() const
  {
    return std::make_unique<ScalarPolynomial<dim>>(*this);
  }

  template class ScalarPolynomial<1>;
  template class ScalarPolynomial<2>;
  template class ScalarPolynomial<3>;

} // namespace Simplex

DEAL_II_NAMESPACE_CLOSE
