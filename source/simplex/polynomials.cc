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

    unsigned int
    compute_n_polynomials_pyramid(const unsigned int dim,
                                  const unsigned int degree)
    {
      if (dim == 3)
        {
          if (degree == 1)
            return 5;
        }

      Assert(false, ExcNotImplemented());

      return 0;
    }

    unsigned int
    compute_n_polynomials_wedge(const unsigned int dim,
                                const unsigned int degree)
    {
      if (dim == 3)
        {
          if (degree == 1)
            return 6;
          if (degree == 2)
            return 18;
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
    (void)p;

    Tensor<2, dim> result;

    if (this->degree() < 2)
      return result;

    if (dim == 1)
      {
        if (i == 0)
          {
            result[0][0] = 4.0;
          }
        else if (i == 1)
          {
            result[0][0] = 4.0;
          }
        else if (i == 2)
          {
            result[0][0] = -8.0;
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
    else if (dim == 2)
      {
        if (i == 0)
          {
            for (unsigned j = 0; j < dim; ++j)
              for (unsigned k = 0; k < dim; ++k)
                result[j][k] = 4.0;
          }
        else if (i == 1)
          {
            result[0][0] = 4.0;
          }
        else if (i == 2)
          {
            result[1][1] = 4.0;
          }
        else if (i == 3)
          {
            result[0][0] = -8.0;
            result[0][1] = result[1][0] = -4.0;
          }
        else if (i == 4)
          {
            result[0][1] = result[1][0] = 4.0;
          }
        else if (i == 5)
          {
            result[1][1] = -8.0;
            result[0][1] = result[1][0] = -4.0;
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
    else if (dim == 3)
      {
        if (i == 0)
          {
            for (unsigned j = 0; j < dim; ++j)
              for (unsigned k = 0; k < dim; ++k)
                result[j][k] = 4.0;
          }
        else if (i == 1)
          {
            result[0][0] = 4.0;
          }
        else if (i == 2)
          {
            result[1][1] = 4.0;
          }
        else if (i == 3)
          {
            result[2][2] = 4.0;
          }
        else if (i == 4)
          {
            result[0][0] = -8.0;
            result[0][1] = result[0][2] = result[1][0] = result[2][0] = -4.0;
          }
        else if (i == 5)
          {
            result[0][1] = result[1][0] = 4.0;
          }
        else if (i == 6)
          {
            result[1][1] = -8.0;
            result[0][1] = result[1][0] = result[1][2] = result[2][1] = -4.0;
          }
        else if (i == 7)
          {
            result[2][2] = -8.0;
            result[2][0] = result[2][1] = result[0][2] = result[1][2] = -4.0;
          }
        else if (i == 8)
          {
            result[0][2] = result[2][0] = 4.0;
          }
        else if (i == 9)
          {
            result[1][2] = result[2][1] = 4.0;
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

    return result;
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

    if (grad_grads.size() == this->n())
      for (unsigned int i = 0; i < this->n(); i++)
        grad_grads[i] = compute_grad_grad(i, unit_point);
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
    return compute_grad_grad(i, p);
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



  template <int dim>
  ScalarWedgePolynomial<dim>::ScalarWedgePolynomial(const unsigned int degree)
    : ScalarPolynomialsBase<dim>(degree,
                                 compute_n_polynomials_wedge(dim, degree))
    , poly_tri(degree)
    , poly_line(degree)
  {}


  namespace
  {
    /**
     * Decompose the shape-function index of a linear wedge into an index
     * to access the right shape function within the triangle and and within
     * the line.
     */
    static const constexpr std::array<std::array<unsigned int, 2>, 6>
      wedge_table_1{
        {{{0, 0}}, {{1, 0}}, {{2, 0}}, {{0, 1}}, {{1, 1}}, {{2, 1}}}};

    /**
     * Decompose the shape-function index of a quadratic wedge into an index
     * to access the right shape function within the triangle and and within
     * the line.
     */
    static const constexpr std::array<std::array<unsigned int, 2>, 18>
      wedge_table_2{{{{0, 0}},
                     {{1, 0}},
                     {{2, 0}},
                     {{0, 1}},
                     {{1, 1}},
                     {{2, 1}},
                     {{3, 0}},
                     {{4, 0}},
                     {{5, 0}},
                     {{3, 1}},
                     {{4, 1}},
                     {{5, 1}},
                     {{0, 2}},
                     {{1, 2}},
                     {{2, 2}},
                     {{3, 2}},
                     {{4, 2}},
                     {{5, 2}}}};
  } // namespace


  template <int dim>
  double
  ScalarWedgePolynomial<dim>::compute_value(const unsigned int i,
                                            const Point<dim> & p) const
  {
    const auto pair = this->degree() == 1 ? wedge_table_1[i] : wedge_table_2[i];

    const Point<2> p_tri(p[0], p[1]);
    const auto     v_tri = poly_tri.compute_value(pair[0], p_tri);

    const Point<1> p_line(p[2]);
    const auto     v_line = poly_line.compute_value(pair[1], p_line);

    return v_tri * v_line;
  }



  template <int dim>
  Tensor<1, dim>
  ScalarWedgePolynomial<dim>::compute_grad(const unsigned int i,
                                           const Point<dim> & p) const
  {
    const auto pair = this->degree() == 1 ? wedge_table_1[i] : wedge_table_2[i];

    const Point<2> p_tri(p[0], p[1]);
    const auto     v_tri = poly_tri.compute_value(pair[0], p_tri);
    const auto     g_tri = poly_tri.compute_grad(pair[0], p_tri);

    const Point<1> p_line(p[2]);
    const auto     v_line = poly_line.compute_value(pair[1], p_line);
    const auto     g_line = poly_line.compute_grad(pair[1], p_line);

    Tensor<1, dim> grad;
    grad[0] = g_tri[0] * v_line;
    grad[1] = g_tri[1] * v_line;
    grad[2] = v_tri * g_line[0];

    return grad;
  }



  template <int dim>
  Tensor<2, dim>
  ScalarWedgePolynomial<dim>::compute_grad_grad(const unsigned int i,
                                                const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());
    return Tensor<2, dim>();
  }



  template <int dim>
  void
  ScalarWedgePolynomial<dim>::evaluate(
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
  ScalarWedgePolynomial<dim>::compute_1st_derivative(const unsigned int i,
                                                     const Point<dim> & p) const
  {
    return compute_grad(i, p);
  }



  template <int dim>
  Tensor<2, dim>
  ScalarWedgePolynomial<dim>::compute_2nd_derivative(const unsigned int i,
                                                     const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());

    return {};
  }



  template <int dim>
  Tensor<3, dim>
  ScalarWedgePolynomial<dim>::compute_3rd_derivative(const unsigned int i,
                                                     const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());

    return {};
  }



  template <int dim>
  Tensor<4, dim>
  ScalarWedgePolynomial<dim>::compute_4th_derivative(const unsigned int i,
                                                     const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());

    return {};
  }



  template <int dim>
  std::string
  ScalarWedgePolynomial<dim>::name() const
  {
    return "ScalarWedgePolynomial";
  }



  template <int dim>
  std::unique_ptr<ScalarPolynomialsBase<dim>>
  ScalarWedgePolynomial<dim>::clone() const
  {
    return std::make_unique<ScalarWedgePolynomial<dim>>(*this);
  }



  template <int dim>
  ScalarPyramidPolynomial<dim>::ScalarPyramidPolynomial(
    const unsigned int degree)
    : ScalarPolynomialsBase<dim>(degree,
                                 compute_n_polynomials_pyramid(dim, degree))
  {}


  template <int dim>
  double
  ScalarPyramidPolynomial<dim>::compute_value(const unsigned int i,
                                              const Point<dim> & p) const
  {
    AssertDimension(dim, 3);
    AssertIndexRange(this->degree(), 2);

    const double Q14 = 0.25;
    double       ration;

    const double r = p[0];
    const double s = p[1];
    const double t = p[2];

    if (fabs(t - 1.0) > 1.0e-14)
      {
        ration = (r * s * t) / (1.0 - t);
      }
    else
      {
        ration = 0.0;
      }

    if (i == 0)
      return Q14 * ((1.0 - r) * (1.0 - s) - t + ration);
    if (i == 1)
      return Q14 * ((1.0 + r) * (1.0 - s) - t - ration);
    if (i == 2)
      return Q14 * ((1.0 - r) * (1.0 + s) - t - ration);
    if (i == 3)
      return Q14 * ((1.0 + r) * (1.0 + s) - t + ration);
    else
      return t;
  }



  template <int dim>
  Tensor<1, dim>
  ScalarPyramidPolynomial<dim>::compute_grad(const unsigned int i,
                                             const Point<dim> & p) const
  {
    AssertDimension(dim, 3);
    AssertIndexRange(this->degree(), 4);

    Tensor<1, dim> grad;

    if (this->degree() == 1)
      {
        const double Q14 = 0.25;

        const double r = p[0];
        const double s = p[1];
        const double t = p[2];

        double rationdr;
        double rationds;
        double rationdt;

        if (fabs(t - 1.0) > 1.0e-14)
          {
            rationdr = s * t / (1.0 - t);
            rationds = r * t / (1.0 - t);
            rationdt = r * s / ((1.0 - t) * (1.0 - t));
          }
        else
          {
            rationdr = 1.0;
            rationds = 1.0;
            rationdt = 1.0;
          }


        if (i == 0)
          {
            grad[0] = Q14 * (-1.0 * (1.0 - s) + rationdr);
            grad[1] = Q14 * (-1.0 * (1.0 - r) + rationds);
            grad[2] = Q14 * (rationdt - 1.0);
          }
        else if (i == 1)
          {
            grad[0] = Q14 * (1.0 * (1.0 - s) - rationdr);
            grad[1] = Q14 * (-1.0 * (1.0 + r) - rationds);
            grad[2] = Q14 * (-1.0 * rationdt - 1.0);
          }
        else if (i == 2)
          {
            grad[0] = Q14 * (-1.0 * (1.0 + s) - rationdr);
            grad[1] = Q14 * (1.0 * (1.0 - r) - rationds);
            grad[2] = Q14 * (-1.0 * rationdt - 1.0);
          }
        else if (i == 3)
          {
            grad[0] = Q14 * (1.0 * (1.0 + s) + rationdr);
            grad[1] = Q14 * (1.0 * (1.0 + r) + rationds);
            grad[2] = Q14 * (rationdt - 1.0);
          }
        else if (i == 4)
          {
            grad[0] = 0.0;
            grad[1] = 0.0;
            grad[2] = 1.0;
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }

    return grad;
  }



  template <int dim>
  Tensor<2, dim>
  ScalarPyramidPolynomial<dim>::compute_grad_grad(const unsigned int i,
                                                  const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());
    return Tensor<2, dim>();
  }



  template <int dim>
  void
  ScalarPyramidPolynomial<dim>::evaluate(
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
  ScalarPyramidPolynomial<dim>::compute_1st_derivative(
    const unsigned int i,
    const Point<dim> & p) const
  {
    return compute_grad(i, p);
  }



  template <int dim>
  Tensor<2, dim>
  ScalarPyramidPolynomial<dim>::compute_2nd_derivative(
    const unsigned int i,
    const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());

    return {};
  }



  template <int dim>
  Tensor<3, dim>
  ScalarPyramidPolynomial<dim>::compute_3rd_derivative(
    const unsigned int i,
    const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());

    return {};
  }



  template <int dim>
  Tensor<4, dim>
  ScalarPyramidPolynomial<dim>::compute_4th_derivative(
    const unsigned int i,
    const Point<dim> & p) const
  {
    (void)i;
    (void)p;

    Assert(false, ExcNotImplemented());

    return {};
  }



  template <int dim>
  std::string
  ScalarPyramidPolynomial<dim>::name() const
  {
    return "ScalarPyramidPolynomial";
  }



  template <int dim>
  std::unique_ptr<ScalarPolynomialsBase<dim>>
  ScalarPyramidPolynomial<dim>::clone() const
  {
    return std::make_unique<ScalarPyramidPolynomial<dim>>(*this);
  }



  template class ScalarPolynomial<1>;
  template class ScalarPolynomial<2>;
  template class ScalarPolynomial<3>;
  template class ScalarWedgePolynomial<1>;
  template class ScalarWedgePolynomial<2>;
  template class ScalarWedgePolynomial<3>;
  template class ScalarPyramidPolynomial<1>;
  template class ScalarPyramidPolynomial<2>;
  template class ScalarPyramidPolynomial<3>;

} // namespace Simplex

DEAL_II_NAMESPACE_CLOSE
