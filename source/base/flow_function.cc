// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/flow_function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <cmath>


DEAL_II_NAMESPACE_OPEN


namespace Functions
{
  template <int dim>
  FlowFunction<dim>::FlowFunction()
    : Function<dim>(dim + 1)
    , mean_pressure(0)
    , aux_values(dim + 1)
    , aux_gradients(dim + 1)
  {}



  template <int dim>
  void
  FlowFunction<dim>::pressure_adjustment(double p)
  {
    mean_pressure = p;
  }


  template <int dim>
  void
  FlowFunction<dim>::vector_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>>   &values) const
  {
    const unsigned int n_points = points.size();
    Assert(values.size() == n_points,
           ExcDimensionMismatch(values.size(), n_points));

    // guard access to the aux_*
    // variables in multithread mode
    std::lock_guard<std::mutex> lock(mutex);

    for (unsigned int d = 0; d < dim + 1; ++d)
      aux_values[d].resize(n_points);
    vector_values(points, aux_values);

    for (unsigned int k = 0; k < n_points; ++k)
      {
        Assert(values[k].size() == dim + 1,
               ExcDimensionMismatch(values[k].size(), dim + 1));
        for (unsigned int d = 0; d < dim + 1; ++d)
          values[k](d) = aux_values[d][k];
      }
  }


  template <int dim>
  void
  FlowFunction<dim>::vector_value(const Point<dim> &point,
                                  Vector<double>   &value) const
  {
    Assert(value.size() == dim + 1,
           ExcDimensionMismatch(value.size(), dim + 1));

    const unsigned int      n_points = 1;
    std::vector<Point<dim>> points(1);
    points[0] = point;

    // guard access to the aux_*
    // variables in multithread mode
    std::lock_guard<std::mutex> lock(mutex);

    for (unsigned int d = 0; d < dim + 1; ++d)
      aux_values[d].resize(n_points);
    vector_values(points, aux_values);

    for (unsigned int d = 0; d < dim + 1; ++d)
      value(d) = aux_values[d][0];
  }


  template <int dim>
  double
  FlowFunction<dim>::value(const Point<dim>  &point,
                           const unsigned int comp) const
  {
    AssertIndexRange(comp, dim + 1);
    const unsigned int      n_points = 1;
    std::vector<Point<dim>> points(1);
    points[0] = point;

    // guard access to the aux_*
    // variables in multithread mode
    std::lock_guard<std::mutex> lock(mutex);

    for (unsigned int d = 0; d < dim + 1; ++d)
      aux_values[d].resize(n_points);
    vector_values(points, aux_values);

    return aux_values[comp][0];
  }


  template <int dim>
  void
  FlowFunction<dim>::vector_gradient_list(
    const std::vector<Point<dim>>            &points,
    std::vector<std::vector<Tensor<1, dim>>> &values) const
  {
    const unsigned int n_points = points.size();
    Assert(values.size() == n_points,
           ExcDimensionMismatch(values.size(), n_points));

    // guard access to the aux_*
    // variables in multithread mode
    std::lock_guard<std::mutex> lock(mutex);

    for (unsigned int d = 0; d < dim + 1; ++d)
      aux_gradients[d].resize(n_points);
    vector_gradients(points, aux_gradients);

    for (unsigned int k = 0; k < n_points; ++k)
      {
        Assert(values[k].size() == dim + 1,
               ExcDimensionMismatch(values[k].size(), dim + 1));
        for (unsigned int d = 0; d < dim + 1; ++d)
          values[k][d] = aux_gradients[d][k];
      }
  }


  template <int dim>
  void
  FlowFunction<dim>::vector_laplacian_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>>   &values) const
  {
    const unsigned int n_points = points.size();
    Assert(values.size() == n_points,
           ExcDimensionMismatch(values.size(), n_points));

    // guard access to the aux_*
    // variables in multithread mode
    std::lock_guard<std::mutex> lock(mutex);

    for (unsigned int d = 0; d < dim + 1; ++d)
      aux_values[d].resize(n_points);
    vector_laplacians(points, aux_values);

    for (unsigned int k = 0; k < n_points; ++k)
      {
        Assert(values[k].size() == dim + 1,
               ExcDimensionMismatch(values[k].size(), dim + 1));
        for (unsigned int d = 0; d < dim + 1; ++d)
          values[k](d) = aux_values[d][k];
      }
  }


  template <int dim>
  std::size_t
  FlowFunction<dim>::memory_consumption() const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }


  //----------------------------------------------------------------------//

  template <int dim>
  PoisseuilleFlow<dim>::PoisseuilleFlow(const double r, const double Re)
    : inv_sqr_radius(1 / r / r)
    , Reynolds(Re)
  {
    Assert(Reynolds != 0., ExcMessage("Reynolds number cannot be zero"));
  }



  template <int dim>
  void
  PoisseuilleFlow<dim>::vector_values(
    const std::vector<Point<dim>>    &points,
    std::vector<std::vector<double>> &values) const
  {
    const unsigned int n = points.size();

    Assert(values.size() == dim + 1,
           ExcDimensionMismatch(values.size(), dim + 1));
    for (unsigned int d = 0; d < dim + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k = 0; k < n; ++k)
      {
        const Point<dim> &p = points[k];
        // First, compute the square of the distance to the x-axis divided by
        // the radius.
        double r2 = 0;
        for (unsigned int d = 1; d < dim; ++d)
          r2 += p[d] * p[d];
        r2 *= inv_sqr_radius;

        // x-velocity
        values[0][k] = 1. - r2;
        // other velocities
        for (unsigned int d = 1; d < dim; ++d)
          values[d][k] = 0.;
        // pressure
        values[dim][k] = -2 * (dim - 1) * inv_sqr_radius * p[0] / Reynolds +
                         this->mean_pressure;
      }
  }



  template <int dim>
  void
  PoisseuilleFlow<dim>::vector_gradients(
    const std::vector<Point<dim>>            &points,
    std::vector<std::vector<Tensor<1, dim>>> &values) const
  {
    const unsigned int n = points.size();

    Assert(values.size() == dim + 1,
           ExcDimensionMismatch(values.size(), dim + 1));
    for (unsigned int d = 0; d < dim + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k = 0; k < n; ++k)
      {
        const Point<dim> &p = points[k];
        // x-velocity
        values[0][k][0] = 0.;
        for (unsigned int d = 1; d < dim; ++d)
          values[0][k][d] = -2. * p[d] * inv_sqr_radius;
        // other velocities
        for (unsigned int d = 1; d < dim; ++d)
          values[d][k] = 0.;
        // pressure
        values[dim][k][0] = -2 * (dim - 1) * inv_sqr_radius / Reynolds;
        for (unsigned int d = 1; d < dim; ++d)
          values[dim][k][d] = 0.;
      }
  }



  template <int dim>
  void
  PoisseuilleFlow<dim>::vector_laplacians(
    const std::vector<Point<dim>>    &points,
    std::vector<std::vector<double>> &values) const
  {
    const unsigned int n = points.size();
    (void)n;
    Assert(values.size() == dim + 1,
           ExcDimensionMismatch(values.size(), dim + 1));
    for (unsigned int d = 0; d < dim + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (auto &point_values : values)
      std::fill(point_values.begin(), point_values.end(), 0.);
  }

  //----------------------------------------------------------------------//

  template <int dim>
  StokesCosine<dim>::StokesCosine(const double nu, const double r)
    : viscosity(nu)
    , reaction(r)
  {}



  template <int dim>
  void
  StokesCosine<dim>::set_parameters(const double nu, const double r)
  {
    viscosity = nu;
    reaction  = r;
  }


  template <int dim>
  void
  StokesCosine<dim>::vector_values(
    const std::vector<Point<dim>>    &points,
    std::vector<std::vector<double>> &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == dim + 1,
           ExcDimensionMismatch(values.size(), dim + 1));
    for (unsigned int d = 0; d < dim + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k = 0; k < n; ++k)
      {
        const Point<dim> &p  = points[k];
        const double      x  = numbers::PI / 2. * p[0];
        const double      y  = numbers::PI / 2. * p[1];
        const double      cx = std::cos(x);
        const double      cy = std::cos(y);
        const double      sx = std::sin(x);
        const double      sy = std::sin(y);

        if (dim == 2)
          {
            values[0][k] = cx * cx * cy * sy;
            values[1][k] = -cx * sx * cy * cy;
            values[2][k] = cx * sx * cy * sy + this->mean_pressure;
          }
        else if (dim == 3)
          {
            const double z  = numbers::PI / 2. * p[2];
            const double cz = std::cos(z);
            const double sz = std::sin(z);

            values[0][k] = cx * cx * cy * sy * cz * sz;
            values[1][k] = cx * sx * cy * cy * cz * sz;
            values[2][k] = -2. * cx * sx * cy * sy * cz * cz;
            values[3][k] = cx * sx * cy * sy * cz * sz + this->mean_pressure;
          }
        else
          {
            DEAL_II_NOT_IMPLEMENTED();
          }
      }
  }



  template <int dim>
  void
  StokesCosine<dim>::vector_gradients(
    const std::vector<Point<dim>>            &points,
    std::vector<std::vector<Tensor<1, dim>>> &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == dim + 1,
           ExcDimensionMismatch(values.size(), dim + 1));
    for (unsigned int d = 0; d < dim + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k = 0; k < n; ++k)
      {
        const Point<dim> &p   = points[k];
        const double      x   = numbers::PI / 2. * p[0];
        const double      y   = numbers::PI / 2. * p[1];
        const double      c2x = std::cos(2 * x);
        const double      c2y = std::cos(2 * y);
        const double      s2x = std::sin(2 * x);
        const double      s2y = std::sin(2 * y);
        const double      cx2 = .5 + .5 * c2x; // cos^2 x
        const double      cy2 = .5 + .5 * c2y; // cos^2 y

        if (dim == 2)
          {
            values[0][k][0] = -.25 * numbers::PI * s2x * s2y;
            values[0][k][1] = .5 * numbers::PI * cx2 * c2y;
            values[1][k][0] = -.5 * numbers::PI * c2x * cy2;
            values[1][k][1] = .25 * numbers::PI * s2x * s2y;
            values[2][k][0] = .25 * numbers::PI * c2x * s2y;
            values[2][k][1] = .25 * numbers::PI * s2x * c2y;
          }
        else if (dim == 3)
          {
            const double z   = numbers::PI / 2. * p[2];
            const double c2z = std::cos(2 * z);
            const double s2z = std::sin(2 * z);
            const double cz2 = .5 + .5 * c2z; // cos^2 z

            values[0][k][0] = -.125 * numbers::PI * s2x * s2y * s2z;
            values[0][k][1] = .25 * numbers::PI * cx2 * c2y * s2z;
            values[0][k][2] = .25 * numbers::PI * cx2 * s2y * c2z;

            values[1][k][0] = .25 * numbers::PI * c2x * cy2 * s2z;
            values[1][k][1] = -.125 * numbers::PI * s2x * s2y * s2z;
            values[1][k][2] = .25 * numbers::PI * s2x * cy2 * c2z;

            values[2][k][0] = -.5 * numbers::PI * c2x * s2y * cz2;
            values[2][k][1] = -.5 * numbers::PI * s2x * c2y * cz2;
            values[2][k][2] = .25 * numbers::PI * s2x * s2y * s2z;

            values[3][k][0] = .125 * numbers::PI * c2x * s2y * s2z;
            values[3][k][1] = .125 * numbers::PI * s2x * c2y * s2z;
            values[3][k][2] = .125 * numbers::PI * s2x * s2y * c2z;
          }
        else
          {
            DEAL_II_NOT_IMPLEMENTED();
          }
      }
  }



  template <int dim>
  void
  StokesCosine<dim>::vector_laplacians(
    const std::vector<Point<dim>>    &points,
    std::vector<std::vector<double>> &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == dim + 1,
           ExcDimensionMismatch(values.size(), dim + 1));
    for (unsigned int d = 0; d < dim + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    if (reaction != 0.)
      {
        vector_values(points, values);
        for (unsigned int d = 0; d < dim; ++d)
          for (double &point_value : values[d])
            point_value *= -reaction;
      }
    else
      {
        for (unsigned int d = 0; d < dim; ++d)
          std::fill(values[d].begin(), values[d].end(), 0.);
      }


    for (unsigned int k = 0; k < n; ++k)
      {
        const Point<dim> &p   = points[k];
        const double      x   = numbers::PI / 2. * p[0];
        const double      y   = numbers::PI / 2. * p[1];
        const double      c2x = std::cos(2 * x);
        const double      c2y = std::cos(2 * y);
        const double      s2x = std::sin(2 * x);
        const double      s2y = std::sin(2 * y);
        const double      pi2 = .25 * numbers::PI * numbers::PI;

        if (dim == 2)
          {
            values[0][k] += -viscosity * pi2 * (1. + 2. * c2x) * s2y -
                            numbers::PI / 4. * c2x * s2y;
            values[1][k] += viscosity * pi2 * s2x * (1. + 2. * c2y) -
                            numbers::PI / 4. * s2x * c2y;
            values[2][k] = 0.;
          }
        else if (dim == 3)
          {
            const double z   = numbers::PI * p[2];
            const double c2z = std::cos(2 * z);
            const double s2z = std::sin(2 * z);

            values[0][k] +=
              -.5 * viscosity * pi2 * (1. + 2. * c2x) * s2y * s2z -
              numbers::PI / 8. * c2x * s2y * s2z;
            values[1][k] += .5 * viscosity * pi2 * s2x * (1. + 2. * c2y) * s2z -
                            numbers::PI / 8. * s2x * c2y * s2z;
            values[2][k] +=
              -.5 * viscosity * pi2 * s2x * s2y * (1. + 2. * c2z) -
              numbers::PI / 8. * s2x * s2y * c2z;
            values[3][k] = 0.;
          }
        else
          {
            DEAL_II_NOT_IMPLEMENTED();
          }
      }
  }


  //----------------------------------------------------------------------//

  const double StokesLSingularity::lambda = 0.54448373678246;

  StokesLSingularity::StokesLSingularity()
    : omega(3. / 2. * numbers::PI)
    , coslo(std::cos(lambda * omega))
    , lp(1. + lambda)
    , lm(1. - lambda)
  {}


  inline double
  StokesLSingularity::Psi(double phi) const
  {
    return coslo * (std::sin(lp * phi) / lp - std::sin(lm * phi) / lm) -
           std::cos(lp * phi) + std::cos(lm * phi);
  }


  inline double
  StokesLSingularity::Psi_1(double phi) const
  {
    return coslo * (std::cos(lp * phi) - std::cos(lm * phi)) +
           lp * std::sin(lp * phi) - lm * std::sin(lm * phi);
  }


  inline double
  StokesLSingularity::Psi_2(double phi) const
  {
    return coslo * (lm * std::sin(lm * phi) - lp * std::sin(lp * phi)) +
           lp * lp * std::cos(lp * phi) - lm * lm * std::cos(lm * phi);
  }


  inline double
  StokesLSingularity::Psi_3(double phi) const
  {
    return coslo *
             (lm * lm * std::cos(lm * phi) - lp * lp * std::cos(lp * phi)) +
           lm * lm * lm * std::sin(lm * phi) -
           lp * lp * lp * std::sin(lp * phi);
  }


  inline double
  StokesLSingularity::Psi_4(double phi) const
  {
    return coslo * (lp * lp * lp * std::sin(lp * phi) -
                    lm * lm * lm * std::sin(lm * phi)) +
           lm * lm * lm * lm * std::cos(lm * phi) -
           lp * lp * lp * lp * std::cos(lp * phi);
  }


  void
  StokesLSingularity::vector_values(
    const std::vector<Point<2>>      &points,
    std::vector<std::vector<double>> &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == 2 + 1, ExcDimensionMismatch(values.size(), 2 + 1));
    for (unsigned int d = 0; d < 2 + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k = 0; k < n; ++k)
      {
        const Point<2> &p = points[k];
        const double    x = p[0];
        const double    y = p[1];

        if ((x < 0) || (y < 0))
          {
            const double phi = std::atan2(y, -x) + numbers::PI;
            const double r2  = x * x + y * y;
            const double rl  = std::pow(r2, lambda / 2.);
            const double rl1 = std::pow(r2, lambda / 2. - .5);
            values[0][k] =
              rl * (lp * std::sin(phi) * Psi(phi) + std::cos(phi) * Psi_1(phi));
            values[1][k] =
              rl * (lp * std::cos(phi) * Psi(phi) - std::sin(phi) * Psi_1(phi));
            values[2][k] = -rl1 * (lp * lp * Psi_1(phi) + Psi_3(phi)) / lm +
                           this->mean_pressure;
          }
        else
          {
            for (unsigned int d = 0; d < 3; ++d)
              values[d][k] = 0.;
          }
      }
  }



  void
  StokesLSingularity::vector_gradients(
    const std::vector<Point<2>>            &points,
    std::vector<std::vector<Tensor<1, 2>>> &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == 2 + 1, ExcDimensionMismatch(values.size(), 2 + 1));
    for (unsigned int d = 0; d < 2 + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k = 0; k < n; ++k)
      {
        const Point<2> &p = points[k];
        const double    x = p[0];
        const double    y = p[1];

        if ((x < 0) || (y < 0))
          {
            const double phi  = std::atan2(y, -x) + numbers::PI;
            const double r2   = x * x + y * y;
            const double r    = std::sqrt(r2);
            const double rl   = std::pow(r2, lambda / 2.);
            const double rl1  = std::pow(r2, lambda / 2. - .5);
            const double rl2  = std::pow(r2, lambda / 2. - 1.);
            const double psi  = Psi(phi);
            const double psi1 = Psi_1(phi);
            const double psi2 = Psi_2(phi);
            const double cosp = std::cos(phi);
            const double sinp = std::sin(phi);

            // Derivatives of u with respect to r, phi
            const double udr = lambda * rl1 * (lp * sinp * psi + cosp * psi1);
            const double udp = rl * (lp * cosp * psi + lp * sinp * psi1 -
                                     sinp * psi1 + cosp * psi2);
            // Derivatives of v with respect to r, phi
            const double vdr = lambda * rl1 * (lp * cosp * psi - sinp * psi1);
            const double vdp = rl * (lp * (cosp * psi1 - sinp * psi) -
                                     cosp * psi1 - sinp * psi2);
            // Derivatives of p with respect to r, phi
            const double pdr =
              -(lambda - 1.) * rl2 * (lp * lp * psi1 + Psi_3(phi)) / lm;
            const double pdp = -rl1 * (lp * lp * psi2 + Psi_4(phi)) / lm;
            values[0][k][0]  = cosp * udr - sinp / r * udp;
            values[0][k][1]  = -sinp * udr - cosp / r * udp;
            values[1][k][0]  = cosp * vdr - sinp / r * vdp;
            values[1][k][1]  = -sinp * vdr - cosp / r * vdp;
            values[2][k][0]  = cosp * pdr - sinp / r * pdp;
            values[2][k][1]  = -sinp * pdr - cosp / r * pdp;
          }
        else
          {
            for (unsigned int d = 0; d < 3; ++d)
              values[d][k] = 0.;
          }
      }
  }



  void
  StokesLSingularity::vector_laplacians(
    const std::vector<Point<2>>      &points,
    std::vector<std::vector<double>> &values) const
  {
    unsigned int n = points.size();
    (void)n;
    Assert(values.size() == 2 + 1, ExcDimensionMismatch(values.size(), 2 + 1));
    for (unsigned int d = 0; d < 2 + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (auto &point_values : values)
      std::fill(point_values.begin(), point_values.end(), 0.);
  }


  //----------------------------------------------------------------------//

  Kovasznay::Kovasznay(double Re, bool stokes)
    : Reynolds(Re)
    , stokes(stokes)
  {
    long double r2 = Reynolds / 2.;
    long double b  = 4 * numbers::PI * numbers::PI;
    long double l  = -b / (r2 + std::sqrt(r2 * r2 + b));
    lbda           = l;
    // mean pressure for a domain
    // spreading from -.5 to 1.5 in
    // x-direction
    p_average = 1 / (8 * l) * (std::exp(3. * l) - std::exp(-l));
  }



  void
  Kovasznay::vector_values(const std::vector<Point<2>>      &points,
                           std::vector<std::vector<double>> &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == 2 + 1, ExcDimensionMismatch(values.size(), 2 + 1));
    for (unsigned int d = 0; d < 2 + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k = 0; k < n; ++k)
      {
        const Point<2> &p   = points[k];
        const double    x   = p[0];
        const double    y   = 2. * numbers::PI * p[1];
        const double    elx = std::exp(lbda * x);

        values[0][k] = 1. - elx * std::cos(y);
        values[1][k] = .5 / numbers::PI * lbda * elx * std::sin(y);
        values[2][k] = -.5 * elx * elx + p_average + this->mean_pressure;
      }
  }


  void
  Kovasznay::vector_gradients(
    const std::vector<Point<2>>            &points,
    std::vector<std::vector<Tensor<1, 2>>> &gradients) const
  {
    unsigned int n = points.size();

    Assert(gradients.size() == 3, ExcDimensionMismatch(gradients.size(), 3));
    Assert(gradients[0].size() == n,
           ExcDimensionMismatch(gradients[0].size(), n));

    for (unsigned int i = 0; i < n; ++i)
      {
        const double x = points[i][0];
        const double y = points[i][1];

        const double elx = std::exp(lbda * x);
        const double cy  = std::cos(2 * numbers::PI * y);
        const double sy  = std::sin(2 * numbers::PI * y);

        // u
        gradients[0][i][0] = -lbda * elx * cy;
        gradients[0][i][1] = 2. * numbers::PI * elx * sy;
        gradients[1][i][0] = lbda * lbda / (2 * numbers::PI) * elx * sy;
        gradients[1][i][1] = lbda * elx * cy;
        // p
        gradients[2][i][0] = -lbda * elx * elx;
        gradients[2][i][1] = 0.;
      }
  }



  void
  Kovasznay::vector_laplacians(const std::vector<Point<2>>      &points,
                               std::vector<std::vector<double>> &values) const
  {
    unsigned int n = points.size();
    Assert(values.size() == 2 + 1, ExcDimensionMismatch(values.size(), 2 + 1));
    for (unsigned int d = 0; d < 2 + 1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    if (stokes)
      {
        const double zp = 2. * numbers::PI;
        for (unsigned int k = 0; k < n; ++k)
          {
            const Point<2> &p   = points[k];
            const double    x   = p[0];
            const double    y   = zp * p[1];
            const double    elx = std::exp(lbda * x);
            const double    u   = 1. - elx * std::cos(y);
            const double    ux  = -lbda * elx * std::cos(y);
            const double    uy  = elx * zp * std::sin(y);
            const double    v   = lbda / zp * elx * std::sin(y);
            const double    vx  = lbda * lbda / zp * elx * std::sin(y);
            const double    vy  = zp * lbda / zp * elx * std::cos(y);

            values[0][k] = u * ux + v * uy;
            values[1][k] = u * vx + v * vy;
            values[2][k] = 0.;
          }
      }
    else
      {
        for (auto &point_values : values)
          std::fill(point_values.begin(), point_values.end(), 0.);
      }
  }

  double
  Kovasznay::lambda() const
  {
    return lbda;
  }



  template class FlowFunction<2>;
  template class FlowFunction<3>;
  template class PoisseuilleFlow<2>;
  template class PoisseuilleFlow<3>;
  template class StokesCosine<2>;
  template class StokesCosine<3>;
} // namespace Functions



DEAL_II_NAMESPACE_CLOSE
