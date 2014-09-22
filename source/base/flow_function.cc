// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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

#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/flow_function.h>
#include <deal.II/lac/vector.h>

#include <cmath>


DEAL_II_NAMESPACE_OPEN


namespace Functions
{

  template<int dim>
  FlowFunction<dim>::FlowFunction()
    :
    Function<dim>(dim+1),
    mean_pressure(0),
    aux_values(dim+1),
    aux_gradients(dim+1)
  {}


  template<int dim>
  FlowFunction<dim>::~FlowFunction()
  {}


  template<int dim>
  void
  FlowFunction<dim>::pressure_adjustment(double p)
  {
    mean_pressure = p;
  }


  template<int dim>
  void FlowFunction<dim>::vector_value_list (
    const std::vector<Point<dim> > &points,
    std::vector<Vector<double> >   &values) const
  {
    const unsigned int n_points = points.size();
    Assert(values.size() == n_points, ExcDimensionMismatch(values.size(), n_points));

    // guard access to the aux_*
    // variables in multithread mode
    Threads::Mutex::ScopedLock lock (mutex);

    for (unsigned int d=0; d<dim+1; ++d)
      aux_values[d].resize(n_points);
    vector_values(points, aux_values);

    for (unsigned int k=0; k<n_points; ++k)
      {
        Assert(values[k].size() == dim+1, ExcDimensionMismatch(values[k].size(), dim+1));
        for (unsigned int d=0; d<dim+1; ++d)
          values[k](d) = aux_values[d][k];
      }
  }


  template<int dim>
  void FlowFunction<dim>::vector_value (
    const Point<dim> &point,
    Vector<double> &value) const
  {
    Assert(value.size() == dim+1, ExcDimensionMismatch(value.size(), dim+1));

    const unsigned int n_points = 1;
    std::vector<Point<dim> > points(1);
    points[0] = point;

    // guard access to the aux_*
    // variables in multithread mode
    Threads::Mutex::ScopedLock lock (mutex);

    for (unsigned int d=0; d<dim+1; ++d)
      aux_values[d].resize(n_points);
    vector_values(points, aux_values);

    for (unsigned int d=0; d<dim+1; ++d)
      value(d) = aux_values[d][0];
  }


  template<int dim>
  double FlowFunction<dim>::value (
    const Point<dim> &point,
    const unsigned int comp) const
  {
    Assert(comp < dim+1, ExcIndexRange(comp, 0, dim+1));
    const unsigned int n_points = 1;
    std::vector<Point<dim> > points(1);
    points[0] = point;

    // guard access to the aux_*
    // variables in multithread mode
    Threads::Mutex::ScopedLock lock (mutex);

    for (unsigned int d=0; d<dim+1; ++d)
      aux_values[d].resize(n_points);
    vector_values(points, aux_values);

    return aux_values[comp][0];
  }


  template<int dim>
  void FlowFunction<dim>::vector_gradient_list (
    const std::vector<Point<dim> > &points,
    std::vector<std::vector<Tensor<1,dim> > > &values) const
  {
    const unsigned int n_points = points.size();
    Assert(values.size() == n_points, ExcDimensionMismatch(values.size(), n_points));

    // guard access to the aux_*
    // variables in multithread mode
    Threads::Mutex::ScopedLock lock (mutex);

    for (unsigned int d=0; d<dim+1; ++d)
      aux_gradients[d].resize(n_points);
    vector_gradients(points, aux_gradients);

    for (unsigned int k=0; k<n_points; ++k)
      {
        Assert(values[k].size() == dim+1, ExcDimensionMismatch(values[k].size(), dim+1));
        for (unsigned int d=0; d<dim+1; ++d)
          values[k][d] = aux_gradients[d][k];
      }
  }


  template<int dim>
  void FlowFunction<dim>::vector_laplacian_list (
    const std::vector<Point<dim> > &points,
    std::vector<Vector<double> > &values) const
  {
    const unsigned int n_points = points.size();
    Assert(values.size() == n_points, ExcDimensionMismatch(values.size(), n_points));

    // guard access to the aux_*
    // variables in multithread mode
    Threads::Mutex::ScopedLock lock (mutex);

    for (unsigned int d=0; d<dim+1; ++d)
      aux_values[d].resize(n_points);
    vector_laplacians(points, aux_values);

    for (unsigned int k=0; k<n_points; ++k)
      {
        Assert(values[k].size() == dim+1, ExcDimensionMismatch(values[k].size(), dim+1));
        for (unsigned int d=0; d<dim+1; ++d)
          values[k](d) = aux_values[d][k];
      }
  }


  template<int dim>
  std::size_t
  FlowFunction<dim>::memory_consumption () const
  {
    Assert(false, ExcNotImplemented());
    return 0;
  }


//----------------------------------------------------------------------//

  template<int dim>
  PoisseuilleFlow<dim>::PoisseuilleFlow(const double r,
                                        const double Re)
    :
    radius(r), Reynolds(Re)
  {
    Assert(Reynolds != 0., ExcMessage("Reynolds number cannot be zero"));
  }


  template<int dim>
  PoisseuilleFlow<dim>::~PoisseuilleFlow()
  {}


  template<int dim>
  void PoisseuilleFlow<dim>::vector_values (
    const std::vector<Point<dim> > &points,
    std::vector<std::vector<double> > &values) const
  {
    unsigned int n = points.size();
    double stretch = 1./radius;

    Assert(values.size() == dim+1, ExcDimensionMismatch(values.size(), dim+1));
    for (unsigned int d=0; d<dim+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k=0; k<n; ++k)
      {
        const Point<dim> &p = points[k];
        // First, compute the
        // square of the distance to
        // the x-axis divided by the
        // radius.
        double r2 = 0;
        for (unsigned int d=1; d<dim; ++d)
          r2 += p(d)*p(d)*stretch*stretch;

        // x-velocity
        values[0][k] = 1.-r2;
        // other velocities
        for (unsigned int d=1; d<dim; ++d)
          values[d][k] = 0.;
        // pressure
        values[dim][k] = -2*(dim-1)*stretch*stretch*p(0)/Reynolds + this->mean_pressure;
      }
  }



  template<int dim>
  void PoisseuilleFlow<dim>::vector_gradients (
    const std::vector<Point<dim> > &points,
    std::vector<std::vector<Tensor<1,dim> > > &values) const
  {
    unsigned int n = points.size();
    double stretch = 1./radius;

    Assert(values.size() == dim+1, ExcDimensionMismatch(values.size(), dim+1));
    for (unsigned int d=0; d<dim+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k=0; k<n; ++k)
      {
        const Point<dim> &p = points[k];
        // x-velocity
        values[0][k][0] = 0.;
        for (unsigned int d=1; d<dim; ++d)
          values[0][k][d] = -2.*p(d)*stretch*stretch;
        // other velocities
        for (unsigned int d=1; d<dim; ++d)
          values[d][k] = 0.;
        // pressure
        values[dim][k][0] = -2*(dim-1)*stretch*stretch/Reynolds;
        for (unsigned int d=1; d<dim; ++d)
          values[dim][k][d] = 0.;
      }
  }



  template<int dim>
  void PoisseuilleFlow<dim>::vector_laplacians (
    const std::vector<Point<dim> > &points,
    std::vector<std::vector<double> > &values) const
  {
    unsigned int n = points.size();
    Assert(values.size() == dim+1, ExcDimensionMismatch(values.size(), dim+1));
    for (unsigned int d=0; d<dim+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int d=0; d<values.size(); ++d)
      for (unsigned int k=0; k<values[d].size(); ++k)
        values[d][k] = 0.;
  }

//----------------------------------------------------------------------//

  template<int dim>
  StokesCosine<dim>::StokesCosine(const double nu, const double r)
    :
    viscosity(nu), reaction(r)
  {}


  template<int dim>
  StokesCosine<dim>::~StokesCosine()
  {}


  template<int dim>
  void
  StokesCosine<dim>::set_parameters(const double nu, const double r)
  {
    viscosity = nu;
    reaction = r;
  }


  template<int dim>
  void StokesCosine<dim>::vector_values (
    const std::vector<Point<dim> > &points,
    std::vector<std::vector<double> > &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == dim+1, ExcDimensionMismatch(values.size(), dim+1));
    for (unsigned int d=0; d<dim+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k=0; k<n; ++k)
      {
        const Point<dim> &p = points[k];
        const double x = numbers::PI/2. * p(0);
        const double y = numbers::PI/2. * p(1);
        const double cx = cos(x);
        const double cy = cos(y);
        const double sx = sin(x);
        const double sy = sin(y);

        if (dim==2)
          {
            values[0][k] = cx*cx*cy*sy;
            values[1][k] = -cx*sx*cy*cy;
            values[2][k] = cx*sx*cy*sy + this->mean_pressure;
          }
        else if (dim==3)
          {
            const double z = numbers::PI/2. * p(2);
            const double cz = cos(z);
            const double sz = sin(z);

            values[0][k] = cx*cx*cy*sy*cz*sz;
            values[1][k] = cx*sx*cy*cy*cz*sz;
            values[2][k] = -2.*cx*sx*cy*sy*cz*cz;
            values[3][k] = cx*sx*cy*sy*cz*sz + this->mean_pressure;
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
  }



  template<int dim>
  void StokesCosine<dim>::vector_gradients (
    const std::vector<Point<dim> > &points,
    std::vector<std::vector<Tensor<1,dim> > > &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == dim+1, ExcDimensionMismatch(values.size(), dim+1));
    for (unsigned int d=0; d<dim+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k=0; k<n; ++k)
      {
        const Point<dim> &p = points[k];
        const double x = numbers::PI/2. * p(0);
        const double y = numbers::PI/2. * p(1);
        const double c2x = cos(2*x);
        const double c2y = cos(2*y);
        const double s2x = sin(2*x);
        const double s2y = sin(2*y);
        const double cx2 = .5+.5*c2x;               // cos^2 x
        const double cy2 = .5+.5*c2y;               // cos^2 y

        if (dim==2)
          {
            values[0][k][0] = -.25*numbers::PI * s2x*s2y;
            values[0][k][1] =  .5 *numbers::PI * cx2*c2y;
            values[1][k][0] = -.5 *numbers::PI * c2x*cy2;
            values[1][k][1] =  .25*numbers::PI * s2x*s2y;
            values[2][k][0] =  .25*numbers::PI * c2x*s2y;
            values[2][k][1] =  .25*numbers::PI * s2x*c2y;
          }
        else if (dim==3)
          {
            const double z = numbers::PI/2. * p(2);
            const double c2z = cos(2*z);
            const double s2z = sin(2*z);
            const double cz2 = .5+.5*c2z;               // cos^2 z

            values[0][k][0] = -.125*numbers::PI * s2x*s2y*s2z;
            values[0][k][1] =  .25 *numbers::PI * cx2*c2y*s2z;
            values[0][k][2] =  .25 *numbers::PI * cx2*s2y*c2z;

            values[1][k][0] =  .25 *numbers::PI * c2x*cy2*s2z;
            values[1][k][1] = -.125*numbers::PI * s2x*s2y*s2z;
            values[1][k][2] =  .25 *numbers::PI * s2x*cy2*c2z;

            values[2][k][0] = -.5  *numbers::PI * c2x*s2y*cz2;
            values[2][k][1] = -.5  *numbers::PI * s2x*c2y*cz2;
            values[2][k][2] =  .25 *numbers::PI * s2x*s2y*s2z;

            values[3][k][0] = .125*numbers::PI * c2x*s2y*s2z;
            values[3][k][1] = .125*numbers::PI * s2x*c2y*s2z;
            values[3][k][2] = .125*numbers::PI * s2x*s2y*c2z;
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
  }



  template<int dim>
  void StokesCosine<dim>::vector_laplacians (
    const std::vector<Point<dim> > &points,
    std::vector<std::vector<double> > &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == dim+1, ExcDimensionMismatch(values.size(), dim+1));
    for (unsigned int d=0; d<dim+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    if (reaction != 0.)
      {
        vector_values(points, values);
        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int k=0; k<values[d].size(); ++k)
            values[d][k] *= -reaction;
      }
    else
      {
        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int k=0; k<values[d].size(); ++k)
            values[d][k] = 0.;
      }


    for (unsigned int k=0; k<n; ++k)
      {
        const Point<dim> &p = points[k];
        const double x = numbers::PI/2. * p(0);
        const double y = numbers::PI/2. * p(1);
        const double c2x = cos(2*x);
        const double c2y = cos(2*y);
        const double s2x = sin(2*x);
        const double s2y = sin(2*y);
        const double pi2 = .25 * numbers::PI * numbers::PI;

        if (dim==2)
          {
            values[0][k] += - viscosity*pi2 * (1.+2.*c2x) * s2y - numbers::PI/4. * c2x*s2y;
            values[1][k] +=   viscosity*pi2 * s2x * (1.+2.*c2y) - numbers::PI/4. * s2x*c2y;
            values[2][k] = 0.;
          }
        else if (dim==3)
          {
            const double z = numbers::PI * p(2);
            const double c2z = cos(2*z);
            const double s2z = sin(2*z);

            values[0][k] += - .5*viscosity*pi2 * (1.+2.*c2x) * s2y * s2z - numbers::PI/8. * c2x * s2y * s2z;
            values[1][k] +=   .5*viscosity*pi2 * s2x * (1.+2.*c2y) * s2z - numbers::PI/8. * s2x * c2y * s2z;
            values[2][k] += - .5*viscosity*pi2 * s2x * s2y * (1.+2.*c2z) - numbers::PI/8. * s2x * s2y * c2z;
            values[3][k] = 0.;
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
  }


//----------------------------------------------------------------------//

  const double StokesLSingularity::lambda = 0.54448373678246;

  StokesLSingularity::StokesLSingularity()
    :
    omega (3./2.*numbers::PI),
    coslo (cos(lambda *omega)),
    lp(1.+lambda),
    lm(1.-lambda)
  {}


  inline
  double
  StokesLSingularity::Psi(double phi) const
  {
    return coslo * (sin(lp*phi)/lp - sin(lm*phi)/lm)
           - cos(lp*phi) + cos(lm*phi);
  }


  inline
  double
  StokesLSingularity::Psi_1(double phi) const
  {
    return coslo * (cos(lp*phi) - cos(lm*phi))
           + lp*sin(lp*phi) - lm*sin(lm*phi);
  }


  inline
  double
  StokesLSingularity::Psi_2(double phi) const
  {
    return coslo * (lm*sin(lm*phi) - lp*sin(lp*phi))
           + lp*lp*cos(lp*phi) - lm*lm*cos(lm*phi);
  }


  inline
  double
  StokesLSingularity::Psi_3(double phi) const
  {
    return coslo * (lm*lm*cos(lm*phi) - lp*lp*cos(lp*phi))
           + lm*lm*lm*sin(lm*phi) - lp*lp*lp*sin(lp*phi);
  }


  inline
  double
  StokesLSingularity::Psi_4(double phi) const
  {
    return coslo * (lp*lp*lp*sin(lp*phi) - lm*lm*lm*sin(lm*phi))
           + lm*lm*lm*lm*cos(lm*phi) - lp*lp*lp*lp*cos(lp*phi);
  }


  void StokesLSingularity::vector_values (
    const std::vector<Point<2> > &points,
    std::vector<std::vector<double> > &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == 2+1, ExcDimensionMismatch(values.size(), 2+1));
    for (unsigned int d=0; d<2+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k=0; k<n; ++k)
      {
        const Point<2> &p = points[k];
        const double x = p(0);
        const double y = p(1);

        if ((x<0) || (y<0))
          {
            const double phi = std::atan2(y,-x)+numbers::PI;
            const double r2 = x*x+y*y;
            const double rl = pow(r2,lambda/2.);
            const double rl1 = pow(r2,lambda/2.-.5);
            values[0][k] = rl * (lp*sin(phi)*Psi(phi) + cos(phi)*Psi_1(phi));
            values[1][k] = rl * (lp*cos(phi)*Psi(phi) - sin(phi)*Psi_1(phi));
            values[2][k] = -rl1 * (lp*lp*Psi_1(phi) + Psi_3(phi)) / lm + this->mean_pressure;
          }
        else
          {
            for (unsigned int d=0; d<3; ++d)
              values[d][k] = 0.;
          }
      }
  }



  void StokesLSingularity::vector_gradients (
    const std::vector<Point<2> > &points,
    std::vector<std::vector<Tensor<1,2> > > &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == 2+1, ExcDimensionMismatch(values.size(), 2+1));
    for (unsigned int d=0; d<2+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k=0; k<n; ++k)
      {
        const Point<2> &p = points[k];
        const double x = p(0);
        const double y = p(1);

        if ((x<0) || (y<0))
          {
            const double phi = std::atan2(y,-x)+numbers::PI;
            const double r2 = x*x+y*y;
            const double r = sqrt(r2);
            const double rl = pow(r2,lambda/2.);
            const double rl1 = pow(r2,lambda/2.-.5);
            const double rl2 = pow(r2,lambda/2.-1.);
            const double psi =Psi(phi);
            const double psi1=Psi_1(phi);
            const double psi2=Psi_2(phi);
            const double cosp= cos(phi);
            const double sinp= sin(phi);

            // Derivatives of u with respect to r, phi
            const double udr = lambda * rl1 * (lp*sinp*psi + cosp*psi1);
            const double udp = rl * (lp*cosp*psi + lp*sinp*psi1 - sinp*psi1 + cosp*psi2);
            // Derivatives of v with respect to r, phi
            const double vdr = lambda * rl1 * (lp*cosp*psi - sinp*psi1);
            const double vdp = rl * (lp*(cosp*psi1 - sinp*psi) - cosp*psi1 - sinp*psi2);
            // Derivatives of p with respect to r, phi
            const double pdr = -(lambda-1.) * rl2 * (lp*lp*psi1+Psi_3(phi)) / lm;
            const double pdp = -rl1 * (lp*lp*psi2+Psi_4(phi)) / lm;
            values[0][k][0] = cosp*udr - sinp/r*udp;
            values[0][k][1] = - sinp*udr - cosp/r*udp;
            values[1][k][0] = cosp*vdr - sinp/r*vdp;
            values[1][k][1] = - sinp*vdr - cosp/r*vdp;
            values[2][k][0] = cosp*pdr - sinp/r*pdp;
            values[2][k][1] = - sinp*pdr - cosp/r*pdp;
          }
        else
          {
            for (unsigned int d=0; d<3; ++d)
              values[d][k] = 0.;
          }
      }
  }



  void StokesLSingularity::vector_laplacians (
    const std::vector<Point<2> > &points,
    std::vector<std::vector<double> > &values) const
  {
    unsigned int n = points.size();
    Assert(values.size() == 2+1, ExcDimensionMismatch(values.size(), 2+1));
    for (unsigned int d=0; d<2+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int d=0; d<values.size(); ++d)
      for (unsigned int k=0; k<values[d].size(); ++k)
        values[d][k] = 0.;
  }


//----------------------------------------------------------------------//

  Kovasznay::Kovasznay(double Re, bool stokes)
    :
    Reynolds(Re),
    stokes(stokes)
  {
    long double r2 = Reynolds/2.;
    long double b = 4*numbers::PI*numbers::PI;
    long double l = -b/(r2+std::sqrt(r2*r2+b));
    lbda = l;
    // mean pressure for a domain
    // spreading from -.5 to 1.5 in
    // x-direction
    p_average = 1/(8*l)*(std::exp(3.*l)-std::exp(-l));
  }


  Kovasznay::~Kovasznay()
  {}


  void Kovasznay::vector_values (
    const std::vector<Point<2> > &points,
    std::vector<std::vector<double> > &values) const
  {
    unsigned int n = points.size();

    Assert(values.size() == 2+1, ExcDimensionMismatch(values.size(), 2+1));
    for (unsigned int d=0; d<2+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    for (unsigned int k=0; k<n; ++k)
      {
        const Point<2> &p = points[k];
        const double x = p(0);
        const double y = 2. * numbers::PI * p(1);
        const double elx = std::exp(lbda*x);

        values[0][k] = 1. - elx * cos(y);
        values[1][k] = .5 / numbers::PI * lbda * elx * sin(y);
        values[2][k] = -.5 * elx * elx + p_average + this->mean_pressure;
      }
  }


  void Kovasznay::vector_gradients (
    const std::vector<Point<2> > &points,
    std::vector<std::vector<Tensor<1,2> > > &gradients) const
  {
    unsigned int n = points.size();

    Assert (gradients.size() == 3, ExcDimensionMismatch(gradients.size(), 3));
    Assert (gradients[0].size() == n,
            ExcDimensionMismatch(gradients[0].size(), n));

    for (unsigned int i=0; i<n; ++i)
      {
        const double x = points[i](0);
        const double y = points[i](1);

        const double elx = std::exp(lbda*x);
        const double cy = cos(2*numbers::PI*y);
        const double sy = sin(2*numbers::PI*y);

        // u
        gradients[0][i][0] = -lbda*elx*cy;
        gradients[0][i][1] = 2. * numbers::PI*elx*sy;
        gradients[1][i][0] = lbda*lbda/(2*numbers::PI)*elx*sy;
        gradients[1][i][1] =lbda*elx*cy;
        // p
        gradients[2][i][0] = -lbda*elx*elx;
        gradients[2][i][1] = 0.;
      }
  }



  void Kovasznay::vector_laplacians (
    const std::vector<Point<2> > &points,
    std::vector<std::vector<double> > &values) const
  {
    unsigned int n = points.size();
    Assert(values.size() == 2+1, ExcDimensionMismatch(values.size(), 2+1));
    for (unsigned int d=0; d<2+1; ++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));

    if (stokes)
      {
        const double zp = 2. * numbers::PI;
        for (unsigned int k=0; k<n; ++k)
          {
            const Point<2> &p = points[k];
            const double x = p(0);
            const double y = zp * p(1);
            const double elx = std::exp(lbda*x);
            const double u  = 1. - elx * cos(y);
            const double ux = -lbda * elx * cos(y);
            const double uy = elx * zp * sin(y);
            const double v  = lbda/zp * elx * sin(y);
            const double vx = lbda*lbda/zp * elx * sin(y);
            const double vy = zp*lbda/zp * elx * cos(y);

            values[0][k] = u*ux+v*uy;
            values[1][k] = u*vx+v*vy;
            values[2][k] = 0.;
          }
      }
    else
      {
        for (unsigned int d=0; d<values.size(); ++d)
          for (unsigned int k=0; k<values[d].size(); ++k)
            values[d][k] = 0.;
      }
  }

  double
  Kovasznay::lambda () const
  {
    return lbda;
  }



  template class FlowFunction<2>;
  template class FlowFunction<3>;
  template class PoisseuilleFlow<2>;
  template class PoisseuilleFlow<3>;
  template class StokesCosine<2>;
  template class StokesCosine<3>;
}



DEAL_II_NAMESPACE_CLOSE
