//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/tensor.h>
#include <base/point.h>
#include <base/flow_function.h>
#include <lac/vector.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


// in strict ANSI C mode, the following constants are not defined by
// default, so we do it ourselves
#ifndef M_PI
#  define	M_PI		3.14159265358979323846
#endif

#ifndef M_PI_2
#  define	M_PI_2		1.57079632679489661923
#endif



namespace Functions
{
  
  template<int dim>
  FlowFunction<dim>::FlowFunction()
		  : Function<dim>(dim+1),
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
    
    for (unsigned int d=0;d<dim+1;++d)
      aux_values[d].resize(n_points);
    vector_values(points, aux_values);

    for (unsigned int k=0;k<n_points;++k)
      {
	Assert(values[k].size() == dim+1, ExcDimensionMismatch(values[k].size(), dim+1));
	for (unsigned int d=0;d<dim+1;++d)
	  values[k](d) = aux_values[d][k];
      }
  }

  
  template<int dim>
  void FlowFunction<dim>::vector_gradient_list (
    const std::vector<Point<dim> >& points,
    std::vector<std::vector<Tensor<1,dim> > >& values) const
  {
    const unsigned int n_points = points.size();
    Assert(values.size() == n_points, ExcDimensionMismatch(values.size(), n_points));
    
    for (unsigned int d=0;d<dim+1;++d)
      aux_gradients[d].resize(n_points);
    vector_gradients(points, aux_gradients);

    for (unsigned int k=0;k<n_points;++k)
      {
	Assert(values[k].size() == dim+1, ExcDimensionMismatch(values[k].size(), dim+1));
	for (unsigned int d=0;d<dim+1;++d)
	  values[k][d] = aux_gradients[d][k];
      }
  }
  
  
  template<int dim>
  void FlowFunction<dim>::vector_laplacian_list (
    const std::vector<Point<dim> >& points,
    std::vector<Vector<double> >& values) const
  {
    const unsigned int n_points = points.size();
    Assert(values.size() == n_points, ExcDimensionMismatch(values.size(), n_points));
    
    for (unsigned int d=0;d<dim+1;++d)
      aux_values[d].resize(n_points);
    vector_laplacians(points, aux_values);

    for (unsigned int k=0;k<n_points;++k)
      {
	Assert(values[k].size() == dim+1, ExcDimensionMismatch(values[k].size(), dim+1));
	for (unsigned int d=0;d<dim+1;++d)
	  values[k](d) = aux_values[d][k];
      }
  }

  
  template<int dim>
  unsigned int FlowFunction<dim>::memory_consumption () const
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
  {}

  
  template<int dim>
  PoisseuilleFlow<dim>::~PoisseuilleFlow()
  {}


  template<int dim>
  void PoisseuilleFlow<dim>::vector_values (
    const std::vector<Point<dim> >& points,
    std::vector<std::vector<double> >& values) const
  {
    unsigned int n = points.size();
    double stretch = 1./radius;

    Assert(values.size() == dim+1, ExcDimensionMismatch(values.size(), dim+1));
    for (unsigned int d=0;d<dim+1;++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));
    
    for (unsigned int k=0;k<n;++k)
      {
	const Point<dim>& p = points[k];
					 // First, compute the
					 // square of the distance to
					 // the x-axis divided by the
					 // radius.
	double r2 = 0;
	for (unsigned int d=1;d<dim;++d)
	  r2 += p(d)*p(d)*stretch*stretch;

					 // x-velocity
	values[0][k] = 1.-r2;
					 // other velocities
	for (unsigned int d=1;d<dim;++d)
	  values[d][k] = 0.;
					 // pressure
	values[dim][k] = -2*(dim-1)*stretch*stretch*p(0)/Reynolds+this->mean_pressure;
      }
  }
  

  
  template<int dim>
  void PoisseuilleFlow<dim>::vector_gradients (
    const std::vector<Point<dim> >& points,
    std::vector<std::vector<Tensor<1,dim> > >& values) const
  {
    unsigned int n = points.size();
    double stretch = 1./radius;

    Assert(values.size() == dim+1, ExcDimensionMismatch(values.size(), dim+1));
    for (unsigned int d=0;d<dim+1;++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));
    
    for (unsigned int k=0;k<n;++k)
      {
	const Point<dim>& p = points[k];
					 // x-velocity
	values[0][k][0] = 0.;
	for (unsigned int d=1;d<dim;++d)
	values[0][k][d] = -2.*p(d)*stretch;
					 // other velocities
	for (unsigned int d=1;d<dim;++d)
	  values[d][k] = 0.;	
					 // pressure
	values[dim][k][0] = -2*(dim-1)*stretch*stretch/Reynolds;
	for (unsigned int d=1;d<dim;++d)
	  values[dim][k][d] = 0.;
      }
  }
  

  
  template<int dim>
  void PoisseuilleFlow<dim>::vector_laplacians (
    const std::vector<Point<dim> >& points,
    std::vector<std::vector<double> >& values) const
  {
    unsigned int n = points.size();
    Assert(values.size() == dim+1, ExcDimensionMismatch(values.size(), dim+1));
    for (unsigned int d=0;d<dim+1;++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));
    
    for (unsigned int d=0;d<values.size();++d)
      for (unsigned int k=0;k<values[d].size();++k)
	values[d][k] = 0.;
  }
  
//----------------------------------------------------------------------//

  const double StokesLSingularity::lambda = 0.54448373678246;
  
  StokesLSingularity::StokesLSingularity()
		  :
		  omega (3./2.*deal_II_numbers::PI),
		  coslo (std::cos(lambda*omega))
  {}
  
  
  inline
  double
  StokesLSingularity::Psi(double phi) const
  {
    return std::sin((1.+lambda) * phi) * coslo / (1.+lambda) - std::cos((1.+lambda) * phi)
      - std::sin((1.-lambda) * phi) * coslo / (1.-lambda) + std::cos((1.-lambda) * phi);
  }
  
  
  inline
  double
  StokesLSingularity::Psi_1(double phi) const
  {
    return std::cos((1.+lambda) * phi) * coslo + (1.+lambda) * std::sin((1.+lambda) * phi)
      - std::cos((1.-lambda) * phi) * coslo - (1.-lambda) * std::sin((1.-lambda) * phi);
  }
  
  
  inline
  double
  StokesLSingularity::Psi_3(double phi) const
  {
    return - (1.+lambda) * (1.+lambda)
      * (std::cos((1.+lambda) * phi) * coslo + (1.+lambda) * std::sin((1.+lambda) * phi))
      - (1.-lambda) * (1.-lambda) *
      (- std::cos((1.-lambda) * phi) * coslo - (1.-lambda) * std::sin((1.-lambda) * phi));
  }
  
  
  void StokesLSingularity::vector_values (
    const std::vector<Point<2> >& points,
    std::vector<std::vector<double> >& values) const
  {
    unsigned int n = points.size();
    
    Assert(values.size() == 2+1, ExcDimensionMismatch(values.size(), 2+1));
    for (unsigned int d=0;d<2+1;++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));
    
    for (unsigned int k=0;k<n;++k)
      {
	const Point<2>& p = points[k];
	const double x = p(0);
	const double y = p(1);

	if ((x<0) || (y<0))
	  {
	    const double phi = std::atan2(y,-x)+M_PI;
	    const double r2 = x*x+y*y;
	    values[0][k] = std::pow(r2,lambda/2.)
			   * ((1.+lambda) * std::sin(phi) * Psi(phi)
			      + std::cos(phi) * Psi_1(phi));
	    values[1][k] = std::pow(r2,lambda/2.)
			   * (std::sin(phi) * Psi_1(phi)
			      -(1.+lambda) * std::cos(phi) * Psi(phi));
	    values[2][k] = -std::pow(r2,lambda/2.-.5)
			   * ((1.+lambda) * (1.+lambda) * Psi_1(phi) + Psi_3(phi))
			   / (1.-lambda);
	  }
	else
	  {
	    for (unsigned int d=0;d<3;++d)
	      values[d][k] = 0.;
	  }
      }
  }
  

  
  void StokesLSingularity::vector_gradients (
    const std::vector<Point<2> >&,
    std::vector<std::vector<Tensor<1,2> > >&) const
  {
    Assert(false, ExcNotImplemented());
  }
  

  
  void StokesLSingularity::vector_laplacians (
    const std::vector<Point<2> >& points,
    std::vector<std::vector<double> >& values) const
  {
    unsigned int n = points.size();
    Assert(values.size() == 2+1, ExcDimensionMismatch(values.size(), 2+1));
    for (unsigned int d=0;d<2+1;++d)
      Assert(values[d].size() == n, ExcDimensionMismatch(values[d].size(), n));
    
    for (unsigned int d=0;d<values.size();++d)
      for (unsigned int k=0;k<values[d].size();++k)
	values[d][k] = 0.;
  }
  
  
  
  template class FlowFunction<2>;
  template class FlowFunction<3>;
  template class PoisseuilleFlow<2>;
  template class PoisseuilleFlow<3>;
}



DEAL_II_NAMESPACE_CLOSE
