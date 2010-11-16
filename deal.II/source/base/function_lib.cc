//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/tensor.h>
#include <base/point.h>
#include <base/function_lib.h>
#include <base/function_bessel.h>
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
  double
  SquareFunction<dim>::value (const Point<dim>   &p,
			      const unsigned int) const
  {
    return p.square();
  }

  
  template<int dim>
  void
  SquareFunction<dim>::vector_value (const Point<dim>   &p,
				     Vector<double>     &values) const
  {
    AssertDimension(values.size(), 1);
    values(0) = p.square();
  }


  template<int dim>
  void
  SquareFunction<dim>::value_list (const std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	values[i] = p.square();
      }
  }
  
  
  template<int dim>
  double
  SquareFunction<dim>::laplacian (const Point<dim>   &,
				  const unsigned int) const
  {
    return 2*dim;
  }
  
  
  template<int dim>
  void
  SquareFunction<dim>::laplacian_list (const std::vector<Point<dim> > &points,
				       std::vector<double>            &values,
				       const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      values[i] = 2*dim;
  }
  
  
  
  template<int dim>
  Tensor<1,dim>
  SquareFunction<dim>::gradient (const Point<dim>   &p,
				 const unsigned int) const
  {
    return p*2.;
  }
  
  
  template<int dim>
  void
  SquareFunction<dim>::vector_gradient (
    const Point<dim>& p,
    std::vector<Tensor<1,dim> >& values) const
  {
    AssertDimension(values.size(), 1);
    values[0] = p*2.;
  }


 
  template<int dim>
  void
  SquareFunction<dim>::gradient_list (const std::vector<Point<dim> > &points,
				      std::vector<Tensor<1,dim> >    &gradients,
				      const unsigned int) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0; i<points.size(); ++i)
      gradients[i] = points[i]*2;
  }
  
  
//////////////////////////////////////////////////////////////////////
  
  
  template<int dim>
  double
  Q1WedgeFunction<dim>::value (const Point<dim>   &p,
			       const unsigned int) const
  {
    Assert (dim>=2, ExcInternalError());
    return p(0)*p(1);
  }
  
  
  
  template<int dim>
  void
  Q1WedgeFunction<dim>::value_list (const std::vector<Point<dim> > &points,
				    std::vector<double>            &values,
				    const unsigned int) const
  {
    Assert (dim>=2, ExcInternalError());
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	values[i] = p(0)*p(1);
      }
  }
  
  
  template<int dim>
  void
  Q1WedgeFunction<dim>::vector_value_list (
    const std::vector<Point<dim> >& points,
    std::vector<Vector<double> >& values) const
  {
    Assert (dim>=2, ExcInternalError());
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    Assert(values[0].size() == 1, ExcDimensionMismatch(values[0].size(), 1));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	values[i](0) = p(0)*p(1);
      }
  }
  
  
  template<int dim>
  double
  Q1WedgeFunction<dim>::laplacian (const Point<dim>   &,
				   const unsigned int) const
  {
    Assert (dim>=2, ExcInternalError());
    return 0.;
  }
  
  
  template<int dim>
  void
  Q1WedgeFunction<dim>::laplacian_list (const std::vector<Point<dim> > &points,
					std::vector<double>            &values,
					const unsigned int) const
  {
    Assert (dim>=2, ExcInternalError());
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      values[i] = 0.;
  }
  
  
  
  template<int dim>
  Tensor<1,dim>
  Q1WedgeFunction<dim>::gradient (const Point<dim>   &p,
				  const unsigned int) const
  {
    Assert (dim>=2, ExcInternalError());
    Tensor<1,dim> erg;
    erg[0] = p(1);
    erg[1] = p(0);
    return erg;
  }
  
  
  
  template<int dim>
  void
  Q1WedgeFunction<dim>::gradient_list (const std::vector<Point<dim> > &points,
				       std::vector<Tensor<1,dim> >    &gradients,
				       const unsigned int) const
  {
    Assert (dim>=2, ExcInternalError());
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0; i<points.size(); ++i)
      {
	gradients[i][0] = points[i](1);
	gradients[i][1] = points[i](0);
      }
  }
  
  
  template<int dim>
  void
  Q1WedgeFunction<dim>::vector_gradient_list (
    const std::vector<Point<dim> >& points,
    std::vector<std::vector<Tensor<1,dim> > >& gradients) const
  {
    Assert (dim>=2, ExcInternalError());
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    Assert(gradients[0].size() == 1,
	   ExcDimensionMismatch(gradients[0].size(), 1));
    
    for (unsigned int i=0; i<points.size(); ++i)
      {
	gradients[i][0][0] = points[i](1);
	gradients[i][0][1] = points[i](0);
      }
  }
  
  
//////////////////////////////////////////////////////////////////////
  
  
  template<int dim>
  PillowFunction<dim>::PillowFunction (const double offset)
		  :
		  offset(offset)
  {}
  
  
  template<int dim>
  double
  PillowFunction<dim>::value (const Point<dim>   &p,
			      const unsigned int) const
  {
    switch(dim)
      {
	case 1:
	      return 1.-p(0)*p(0)+offset;
	case 2:
	      return (1.-p(0)*p(0))*(1.-p(1)*p(1))+offset;
	case 3:
	      return (1.-p(0)*p(0))*(1.-p(1)*p(1))*(1.-p(2)*p(2))+offset;
	default:
	      Assert(false, ExcNotImplemented());
      }
    return 0.;
  }
  
  template<int dim>
  void
  PillowFunction<dim>::value_list (const std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  values[i] = 1.-p(0)*p(0)+offset;
		  break;
	    case 2:
		  values[i] = (1.-p(0)*p(0))*(1.-p(1)*p(1))+offset;
		  break;
	    case 3:
		  values[i] = (1.-p(0)*p(0))*(1.-p(1)*p(1))*(1.-p(2)*p(2))+offset;
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
  
  
  
  template<int dim>
  double
  PillowFunction<dim>::laplacian (const Point<dim>   &p,
				  const unsigned int) const
  {
    switch(dim)
      {
	case 1:
	      return -2.;
	case 2:
	      return -2.*((1.-p(0)*p(0))+(1.-p(1)*p(1)));
	case 3:
	      return -2.*((1.-p(0)*p(0))*(1.-p(1)*p(1))
			  +(1.-p(1)*p(1))*(1.-p(2)*p(2))
			  +(1.-p(2)*p(2))*(1.-p(0)*p(0)));
	default:
	      Assert(false, ExcNotImplemented());
      }
    return 0.;
  }
  
  template<int dim>
  void
  PillowFunction<dim>::laplacian_list (const std::vector<Point<dim> > &points,
				       std::vector<double>            &values,
				       const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  values[i] = -2.;
		  break;
	    case 2:
		  values[i] = -2.*((1.-p(0)*p(0))+(1.-p(1)*p(1)));
		  break;
	    case 3:
		  values[i] = -2.*((1.-p(0)*p(0))*(1.-p(1)*p(1))
				   +(1.-p(1)*p(1))*(1.-p(2)*p(2))
				   +(1.-p(2)*p(2))*(1.-p(0)*p(0)));
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
  
  template<int dim>
  Tensor<1,dim>
  PillowFunction<dim>::gradient (const Point<dim>   &p,
				 const unsigned int) const
  {
    Tensor<1,dim> result;
    switch(dim)
      {
	case 1:
	      result[0] = -2.*p(0);
	      break;
	case 2:
	      result[0] = -2.*p(0)*(1.-p(1)*p(1));
	      result[1] = -2.*p(1)*(1.-p(0)*p(0));
	      break;
	case 3:
	      result[0] = -2.*p(0)*(1.-p(1)*p(1))*(1.-p(2)*p(2));
	      result[1] = -2.*p(1)*(1.-p(0)*p(0))*(1.-p(2)*p(2));
	      result[2] = -2.*p(2)*(1.-p(0)*p(0))*(1.-p(1)*p(1));
	      break;
	default:
	      Assert(false, ExcNotImplemented());
      }
    return result;
  }
  
  template<int dim>
  void
  PillowFunction<dim>::gradient_list (const std::vector<Point<dim> > &points,
				      std::vector<Tensor<1,dim> >    &gradients,
				      const unsigned int) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  gradients[i][0] = -2.*p(0);
		  break;
	    case 2:
		  gradients[i][0] = -2.*p(0)*(1.-p(1)*p(1));
		  gradients[i][1] = -2.*p(1)*(1.-p(0)*p(0));
		  break;
	    case 3:
		  gradients[i][0] = -2.*p(0)*(1.-p(1)*p(1))*(1.-p(2)*p(2));
		  gradients[i][1] = -2.*p(1)*(1.-p(0)*p(0))*(1.-p(2)*p(2));
		  gradients[i][2] = -2.*p(2)*(1.-p(0)*p(0))*(1.-p(1)*p(1));
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
  
//////////////////////////////////////////////////////////////////////

  template <int dim>
  CosineFunction<dim>::CosineFunction (const unsigned int n_components)
		  :
		  Function<dim> (n_components)
  {}



  template<int dim>
  double
  CosineFunction<dim>::value (const Point<dim>   &p,
			      const unsigned int) const
  {
    switch(dim)
      {
	case 1:
	      return std::cos(M_PI_2*p(0));
	case 2:
	      return std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
	case 3:
	      return std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
	default:
	      Assert(false, ExcNotImplemented());
      }
    return 0.;
  }
  
  template<int dim>
  void
  CosineFunction<dim>::value_list (const std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      values[i] = value(points[i]);
  }
  
  
  template<int dim>
  void
  CosineFunction<dim>::vector_value_list (
    const std::vector<Point<dim> > &points,
    std::vector<Vector<double> >   &values) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const double v = value(points[i]);
	for (unsigned int k=0;k<values[i].size();++k)
	  values[i](k) = v;
      }
  }
  
  
  template<int dim>
  double
  CosineFunction<dim>::laplacian (const Point<dim>   &p,
				  const unsigned int) const
  {
    switch(dim)
      {
	case 1:
	      return -M_PI_2*M_PI_2* std::cos(M_PI_2*p(0));
	case 2:
	      return -2*M_PI_2*M_PI_2* std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
	case 3:
	      return -3*M_PI_2*M_PI_2* std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
	default:
	      Assert(false, ExcNotImplemented());
      }
    return 0.;
  }
  
  template<int dim>
  void
  CosineFunction<dim>::laplacian_list (const std::vector<Point<dim> > &points,
				       std::vector<double>            &values,
				       const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      values[i] = laplacian(points[i]);
  }
  
  template<int dim>
  Tensor<1,dim>
  CosineFunction<dim>::gradient (const Point<dim>   &p,
				 const unsigned int) const
  {
    Tensor<1,dim> result;
    switch(dim)
      {
	case 1:
	      result[0] = -M_PI_2* std::sin(M_PI_2*p(0));
	      break;
	case 2:
	      result[0] = -M_PI_2* std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
	      result[1] = -M_PI_2* std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1));
	      break;
	case 3:
	      result[0] = -M_PI_2* std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
	      result[1] = -M_PI_2* std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
	      result[2] = -M_PI_2* std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::sin(M_PI_2*p(2));
	      break;
	default:
	      Assert(false, ExcNotImplemented());
      }
    return result;
  }
  
  template<int dim>
  void
  CosineFunction<dim>::gradient_list (const std::vector<Point<dim> > &points,
				      std::vector<Tensor<1,dim> >    &gradients,
				      const unsigned int) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  gradients[i][0] = -M_PI_2* std::sin(M_PI_2*p(0));
		  break;
	    case 2:
		  gradients[i][0] = -M_PI_2* std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
		  gradients[i][1] = -M_PI_2* std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1));
		  break;
	    case 3:
		  gradients[i][0] = -M_PI_2* std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		  gradients[i][1] = -M_PI_2* std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		  gradients[i][2] = -M_PI_2* std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::sin(M_PI_2*p(2));
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
  
  template<int dim>
  Tensor<2,dim>
  CosineFunction<dim>::hessian (const Point<dim>   &p,
				const unsigned int) const
  {
    const double pi2 = M_PI_2*M_PI_2;
    
    Tensor<2,dim> result;
    switch(dim)
      {
	case 1:
	      result[0][0] = -pi2* std::cos(M_PI_2*p(0));
	      break;
	case 2:
	      if (true)
		{
		  const double coco = -pi2*std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
		  const double sisi = pi2*std::sin(M_PI_2*p(0)) * std::sin(M_PI_2*p(1));
		  result[0][0] = coco;
		  result[1][1] = coco;
		  result[0][1] = sisi;
		  result[1][0] = sisi;
		}
	      break;
	case 3:
	      if (true)
		{
		  const double cococo = -pi2*std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		  const double sisico = pi2*std::sin(M_PI_2*p(0)) * std::sin(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		  const double sicosi = pi2*std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::sin(M_PI_2*p(2));
		  const double cosisi = pi2*std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1)) * std::sin(M_PI_2*p(2));
		  
		  result[0][0] = cococo;
		  result[1][1] = cococo;
		  result[2][2] = cococo;
		  result[0][1] = sisico;
		  result[1][0] = sisico;
		  result[0][2] = sicosi;
		  result[2][0] = sicosi;
		  result[1][2] = cosisi;
		  result[2][1] = cosisi;
		}
	      break;
	default:
	      Assert(false, ExcNotImplemented());
      }
    return result;
  }
  
  template<int dim>
  void
  CosineFunction<dim>::hessian_list (const std::vector<Point<dim> > &points,
				     std::vector<Tensor<2,dim> >    &hessians,
				     const unsigned int) const
  {
    Assert (hessians.size() == points.size(),
	    ExcDimensionMismatch(hessians.size(), points.size()));
    
    const double pi2 = M_PI_2*M_PI_2;
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  hessians[i][0][0] = -pi2* std::cos(M_PI_2*p(0));
		  break;
	    case 2:
		  if (true)
		    {
		      const double coco = -pi2*std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
		      const double sisi = pi2*std::sin(M_PI_2*p(0)) * std::sin(M_PI_2*p(1));
		      hessians[i][0][0] = coco;
		      hessians[i][1][1] = coco;
		      hessians[i][0][1] = sisi;
		      hessians[i][1][0] = sisi;
		    }
		  break;
	    case 3:
		  if (true)
		    {
		      const double cococo = -pi2*std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		      const double sisico = pi2*std::sin(M_PI_2*p(0)) * std::sin(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		      const double sicosi = pi2*std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::sin(M_PI_2*p(2));
		      const double cosisi = pi2*std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1)) * std::sin(M_PI_2*p(2));
		      
		      hessians[i][0][0] = cococo;
		      hessians[i][1][1] = cococo;
		      hessians[i][2][2] = cococo;
		      hessians[i][0][1] = sisico;
		      hessians[i][1][0] = sisico;
		      hessians[i][0][2] = sicosi;
		      hessians[i][2][0] = sicosi;
		      hessians[i][1][2] = cosisi;
		      hessians[i][2][1] = cosisi;
		    }
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
  
//////////////////////////////////////////////////////////////////////
  
  template <int dim>
  CosineGradFunction<dim>::CosineGradFunction ()
		  :
		  Function<dim> (dim)
  {}
  
  
  template<int dim>
  double
  CosineGradFunction<dim>::value (
    const Point<dim>   &p,
    const unsigned int d) const
  {
    AssertIndexRange(d, dim);
    const unsigned int d1 = (d+1) % dim;
    const unsigned int d2 = (d+2) % dim;
    switch(dim)
      {
	case 1:
	      return (-M_PI_2* std::sin(M_PI_2*p(0)));
	case 2:
	      return (-M_PI_2* std::sin(M_PI_2*p(d)) * std::cos(M_PI_2*p(d1)));
	case 3:
	      return (-M_PI_2* std::sin(M_PI_2*p(d)) * std::cos(M_PI_2*p(d1)) * std::cos(M_PI_2*p(d2)));
	default:
	      Assert(false, ExcNotImplemented());
      }
    return 0.;
  }

  
  template<int dim>
  void
  CosineGradFunction<dim>::vector_value (
    const Point<dim>& p,
    Vector<double>& result) const
  {
    AssertDimension(result.size(), dim);
    switch(dim)
      {
	case 1:
	      result(0) = -M_PI_2* std::sin(M_PI_2*p(0));
	      break;
	case 2:
	      result(0) = -M_PI_2* std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
	      result(1) = -M_PI_2* std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1));
	      break;
	case 3:
	      result(0) = -M_PI_2* std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
	      result(1) = -M_PI_2* std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
	      result(2) = -M_PI_2* std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::sin(M_PI_2*p(2));
	      break;
	default:
	      Assert(false, ExcNotImplemented());
      }
  }
  

  template<int dim>
  void
  CosineGradFunction<dim>::value_list (
    const std::vector<Point<dim> >& points,
    std::vector<double>& values,
    const unsigned int d) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    AssertIndexRange(d, dim);
    const unsigned int d1 = (d+1) % dim;
    const unsigned int d2 = (d+2) % dim;
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  values[i] = -M_PI_2* std::sin(M_PI_2*p(d));
		  break;
	    case 2:
		  values[i] = -M_PI_2* std::sin(M_PI_2*p(d)) * std::cos(M_PI_2*p(d1));
		  break;
	    case 3:
		  values[i] = -M_PI_2* std::sin(M_PI_2*p(d)) * std::cos(M_PI_2*p(d1)) * std::cos(M_PI_2*p(d2));
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
  
  
  template<int dim>
  void
  CosineGradFunction<dim>::vector_value_list (
    const std::vector<Point<dim> > &points,
    std::vector<Vector<double> >   &values) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  values[i](0) = -M_PI_2* std::sin(M_PI_2*p(0));
		  break;
	    case 2:
		  values[i](0) = -M_PI_2* std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
		  values[i](1) = -M_PI_2* std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1));
		  break;
	    case 3:
		  values[i](0) = -M_PI_2* std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		  values[i](1) = -M_PI_2* std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		  values[i](2) = -M_PI_2* std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::sin(M_PI_2*p(2));
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
  
  
  template<int dim>
  double
  CosineGradFunction<dim>::laplacian (
    const Point<dim>   &p,
    const unsigned int d) const
  {
    return -M_PI_2*M_PI_2* value(p,d);
  }


  template<int dim>
  Tensor<1,dim>
  CosineGradFunction<dim>::gradient (
    const Point<dim>& p,
    const unsigned int d) const
  {
    AssertIndexRange(d, dim);
    const unsigned int d1 = (d+1) % dim;
    const unsigned int d2 = (d+2) % dim;
    const double pi2 = M_PI_2*M_PI_2;
    
    Tensor<1,dim> result;
    switch(dim)
      {
	case 1:
	      result[0] = -pi2* std::cos(M_PI_2*p(0));
	      break;
	case 2:
	      result[d ] = -pi2*std::cos(M_PI_2*p(d)) * std::cos(M_PI_2*p(d1));
	      result[d1] =  pi2*std::sin(M_PI_2*p(d)) * std::sin(M_PI_2*p(d1));
	      break;
	case 3:
	      result[d ] = -pi2*std::cos(M_PI_2*p(d)) * std::cos(M_PI_2*p(d1)) * std::cos(M_PI_2*p(d2));
	      result[d1] =  pi2*std::sin(M_PI_2*p(d)) * std::sin(M_PI_2*p(d1)) * std::cos(M_PI_2*p(d2));
	      result[d2] =  pi2*std::sin(M_PI_2*p(d)) * std::cos(M_PI_2*p(d1)) * std::sin(M_PI_2*p(d2));
	      break;
	default:
	      Assert(false, ExcNotImplemented());
      }
    return result;
  }

  
  template<int dim>
  void
  CosineGradFunction<dim>::gradient_list (
    const std::vector<Point<dim> > &points,
    std::vector<Tensor<1,dim> >    &gradients,
    const unsigned int d) const
  {
    AssertIndexRange(d, dim);
    const unsigned int d1 = (d+1) % dim;
    const unsigned int d2 = (d+2) % dim;
    const double pi2 = M_PI_2*M_PI_2;
    
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
     for (unsigned int i=0;i<points.size();++i)
       {
	 const Point<dim>& p = points[i];
	 Tensor<1,dim>& result = gradients[i];
	 
	 switch(dim)
	   {
	     case 1:
		   result[0] = -pi2* std::cos(M_PI_2*p(0));
		   break;
	     case 2:
		   result[d ] = -pi2*std::cos(M_PI_2*p(d)) * std::cos(M_PI_2*p(d1));
		   result[d1] =  pi2*std::sin(M_PI_2*p(d)) * std::sin(M_PI_2*p(d1));
		   break;
	     case 3:
		   result[d ] = -pi2*std::cos(M_PI_2*p(d)) * std::cos(M_PI_2*p(d1)) * std::cos(M_PI_2*p(d2));
		   result[d1] =  pi2*std::sin(M_PI_2*p(d)) * std::sin(M_PI_2*p(d1)) * std::cos(M_PI_2*p(d2));
		   result[d2] =  pi2*std::sin(M_PI_2*p(d)) * std::cos(M_PI_2*p(d1)) * std::sin(M_PI_2*p(d2));
		   break;
	     default:
		   Assert(false, ExcNotImplemented());
	   }
       }
  }
  
  
  template<int dim>
  void
  CosineGradFunction<dim>::vector_gradient_list (
    const std::vector<Point<dim> > &points,
    std::vector<std::vector<Tensor<1,dim> > >& gradients) const
  {
    AssertVectorVectorDimension(gradients, points.size(), dim);      
    const double pi2 = M_PI_2*M_PI_2;
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  gradients[i][0][0] = -pi2* std::cos(M_PI_2*p(0));
		  break;
	    case 2:
		  if (true)
		    {
		      const double coco = -pi2*std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
		      const double sisi = pi2*std::sin(M_PI_2*p(0)) * std::sin(M_PI_2*p(1));
		      gradients[i][0][0] = coco;
		      gradients[i][1][1] = coco;
		      gradients[i][0][1] = sisi;
		      gradients[i][1][0] = sisi;
		    }
		  break;
	    case 3:
		  if (true)
		    {
		      const double cococo = -pi2*std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		      const double sisico = pi2*std::sin(M_PI_2*p(0)) * std::sin(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		      const double sicosi = pi2*std::sin(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::sin(M_PI_2*p(2));
		      const double cosisi = pi2*std::cos(M_PI_2*p(0)) * std::sin(M_PI_2*p(1)) * std::sin(M_PI_2*p(2));
			
		      gradients[i][0][0] = cococo;
		      gradients[i][1][1] = cococo;
		      gradients[i][2][2] = cococo;
		      gradients[i][0][1] = sisico;
		      gradients[i][1][0] = sisico;
		      gradients[i][0][2] = sicosi;
		      gradients[i][2][0] = sicosi;
		      gradients[i][1][2] = cosisi;
		      gradients[i][2][1] = cosisi;
		    }
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
    
    
//////////////////////////////////////////////////////////////////////
  
  template<int dim>
  double
  ExpFunction<dim>::value (const Point<dim>   &p,
			   const unsigned int) const
  {
    switch(dim)
      {
	case 1:
	      return std::exp(p(0));
	case 2:
	      return std::exp(p(0)) * std::exp(p(1));
	case 3:
	      return std::exp(p(0)) * std::exp(p(1)) * std::exp(p(2));
	default:
	      Assert(false, ExcNotImplemented());
      }
    return 0.;
  }
  
  template<int dim>
  void
  ExpFunction<dim>::value_list (const std::vector<Point<dim> > &points,
				std::vector<double>            &values,
				const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  values[i] = std::exp(p(0));
		  break;
	    case 2:
		  values[i] = std::exp(p(0)) * std::exp(p(1));
		  break;
	    case 3:
		  values[i] = std::exp(p(0)) * std::exp(p(1)) * std::exp(p(2));
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
  
  template<int dim>
  double
  ExpFunction<dim>::laplacian (const Point<dim>   &p,
			       const unsigned int) const
  {
    switch(dim)
      {
	case 1:
	      return std::exp(p(0));
	case 2:
	      return 2 * std::exp(p(0)) * std::exp(p(1));
	case 3:
	      return 3 * std::exp(p(0)) * std::exp(p(1)) * std::exp(p(2));
	default:
	      Assert(false, ExcNotImplemented());
      }
    return 0.;
  }
  
  template<int dim>
  void
  ExpFunction<dim>::laplacian_list (const std::vector<Point<dim> > &points,
				    std::vector<double>            &values,
				    const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  values[i] = std::exp(p(0));
		  break;
	    case 2:
		  values[i] = 2 * std::exp(p(0)) * std::exp(p(1));
		  break;
	    case 3:
		  values[i] = 3 * std::exp(p(0)) * std::exp(p(1)) * std::exp(p(2));
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
  
  template<int dim>
  Tensor<1,dim>
  ExpFunction<dim>::gradient (const Point<dim>   &p,
			      const unsigned int) const
  {
    Tensor<1,dim> result;
    switch(dim)
      {
	case 1:
	      result[0] = std::exp(p(0));
	      break;
	case 2:
	      result[0] = std::exp(p(0)) * std::exp(p(1));
	      result[1] = std::exp(p(0)) * std::exp(p(1));
	      break;
	case 3:
	      result[0] = std::exp(p(0)) * std::exp(p(1)) * std::exp(p(2));
	      result[1] = std::exp(p(0)) * std::exp(p(1)) * std::exp(p(2));
	      result[2] = std::exp(p(0)) * std::exp(p(1)) * std::exp(p(2));
	      break;
	default:
	      Assert(false, ExcNotImplemented());
      }
    return result;
  }
  
  template<int dim>
  void
  ExpFunction<dim>::gradient_list (const std::vector<Point<dim> > &points,
				   std::vector<Tensor<1,dim> >    &gradients,
				   const unsigned int) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	switch(dim)
	  {
	    case 1:
		  gradients[i][0] = std::exp(p(0));
		  break;
	    case 2:
		  gradients[i][0] = std::exp(p(0)) * std::exp(p(1));
		  gradients[i][1] = std::exp(p(0)) * std::exp(p(1));
		  break;
	    case 3:
		  gradients[i][0] = std::exp(p(0)) * std::exp(p(1)) * std::exp(p(2));
		  gradients[i][1] = std::exp(p(0)) * std::exp(p(1)) * std::exp(p(2));
		  gradients[i][2] = std::exp(p(0)) * std::exp(p(1)) * std::exp(p(2));
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
  }
  
//////////////////////////////////////////////////////////////////////
  
  
  double
  LSingularityFunction::value (const Point<2>   &p,
			       const unsigned int) const
  {
    double x = p(0);
    double y = p(1);
    
    if ((x>=0) && (y>=0))
      return 0.;
    
    double phi = std::atan2(y,-x)+M_PI;
    double r2 = x*x+y*y;
    
    return std::pow(r2,1./3.) * std::sin(2./3.*phi);
  }
  
  
  void
  LSingularityFunction::value_list (const std::vector<Point<2> > &points,
				    std::vector<double>            &values,
				    const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	double x = points[i](0);
	double y = points[i](1);
	
	if ((x>=0) && (y>=0))
	  values[i] = 0.;
	else
	  {
	    double phi = std::atan2(y,-x)+M_PI;
	    double r2 = x*x+y*y;
	    
	    values[i] = std::pow(r2,1./3.) * std::sin(2./3.*phi);
	  }
      }
  }
  
  
  void
  LSingularityFunction::vector_value_list (
    const std::vector<Point<2> >& points,
    std::vector<Vector<double> >& values) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	Assert (values[i].size() == 1,
		ExcDimensionMismatch(values[i].size(), 1));
	double x = points[i](0);
	double y = points[i](1);
	
	if ((x>=0) && (y>=0))
	  values[i](0) = 0.;
	else
	  {
	    double phi = std::atan2(y,-x)+M_PI;
	    double r2 = x*x+y*y;
	    
	    values[i](0) = std::pow(r2,1./3.) * std::sin(2./3.*phi);
	  }
      }
  }
  
  
  double
  LSingularityFunction::laplacian (const Point<2>   &,
				   const unsigned int) const
  {
    return 0.;
  }
  
  
  void
  LSingularityFunction::laplacian_list (const std::vector<Point<2> > &points,
					std::vector<double>            &values,
					const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      values[i] = 0.;
  }
  

  Tensor<1,2>
  LSingularityFunction::gradient (const Point<2>   &p,
				  const unsigned int) const
  {
    double x = p(0);
    double y = p(1);
    double phi = std::atan2(y,-x)+M_PI;
    double r43 = std::pow(x*x+y*y,2./3.);
    
    Tensor<1,2> result;
    result[0] = 2./3.*(std::sin(2./3.*phi)*x + std::cos(2./3.*phi)*y)/r43;
    result[1] = 2./3.*(std::sin(2./3.*phi)*y - std::cos(2./3.*phi)*x)/r43;
    return result;
  }
  
  
  void
  LSingularityFunction::gradient_list (const std::vector<Point<2> > &points,
				       std::vector<Tensor<1,2> >    &gradients,
				       const unsigned int) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<2>& p = points[i];
	double x = p(0);
	double y = p(1);
	double phi = std::atan2(y,-x)+M_PI;
	double r43 = std::pow(x*x+y*y,2./3.);
	
	gradients[i][0] = 2./3.*(std::sin(2./3.*phi)*x + std::cos(2./3.*phi)*y)/r43;
	gradients[i][1] = 2./3.*(std::sin(2./3.*phi)*y - std::cos(2./3.*phi)*x)/r43;
      }
  }
  
  
  void
  LSingularityFunction::vector_gradient_list (
    const std::vector<Point<2> >& points,
    std::vector<std::vector<Tensor<1,2> > >& gradients) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	Assert(gradients[i].size() == 1,
	       ExcDimensionMismatch(gradients[i].size(), 1));
	const Point<2>& p = points[i];
	double x = p(0);
	double y = p(1);
	double phi = std::atan2(y,-x)+M_PI;
	double r43 = std::pow(x*x+y*y,2./3.);
	
	gradients[i][0][0] = 2./3.*(std::sin(2./3.*phi)*x + std::cos(2./3.*phi)*y)/r43;
	gradients[i][0][1] = 2./3.*(std::sin(2./3.*phi)*y - std::cos(2./3.*phi)*x)/r43;
      }
  }
  
//////////////////////////////////////////////////////////////////////
  
  LSingularityGradFunction::LSingularityGradFunction ()
		  :
		  Function<2> (2)
  {}

  
  double
  LSingularityGradFunction::value (const Point<2>   &p,
				   const unsigned int d) const
  {
    AssertIndexRange(d,2);
    
    const double x = p(0);
    const double y = p(1);
    const double phi = std::atan2(y,-x)+M_PI;
    const double r43 = std::pow(x*x+y*y,2./3.);

    return 2./3.*(std::sin(2./3.*phi)*p(d) +
		  (d==0
		   ? (std::cos(2./3.*phi)*y)
		   : (-std::cos(2./3.*phi)*x)))
      /r43;
  }
  
  
  void
  LSingularityGradFunction::value_list (
    const std::vector<Point<2> >& points,
    std::vector<double>& values,
    const unsigned int d) const
  {
    AssertIndexRange(d, 2);
    AssertDimension(values.size(), points.size());
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<2>& p = points[i];
	const double x = p(0);
	const double y = p(1);
	const double phi = std::atan2(y,-x)+M_PI;
	const double r43 = std::pow(x*x+y*y,2./3.);
	
	values[i] = 2./3.*(std::sin(2./3.*phi)*p(d) +
			   (d==0
			    ? (std::cos(2./3.*phi)*y)
			    : (-std::cos(2./3.*phi)*x)))
		    /r43;
      }
  }
  
  
  void
  LSingularityGradFunction::vector_value_list (
    const std::vector<Point<2> >& points,
    std::vector<Vector<double> >& values) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	AssertDimension(values[i].size(), 2);
	const Point<2>& p = points[i];
	const double x = p(0);
	const double y = p(1);
	const double phi = std::atan2(y,-x)+M_PI;
	const double r43 = std::pow(x*x+y*y,2./3.);

	values[i](0) = 2./3.*(std::sin(2./3.*phi)*x + std::cos(2./3.*phi)*y)/r43;
	values[i](1) = 2./3.*(std::sin(2./3.*phi)*y - std::cos(2./3.*phi)*x)/r43;
      }
  }
  
  
  double
  LSingularityGradFunction::laplacian (const Point<2>   &,
				       const unsigned int) const
  {
    return 0.;
  }
  
  
  void
  LSingularityGradFunction::laplacian_list (const std::vector<Point<2> > &points,
					    std::vector<double>            &values,
					    const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      values[i] = 0.;
  }
  
  
  
  Tensor<1,2>
  LSingularityGradFunction::gradient (
    const Point<2>   &/*p*/,
    const unsigned int /*component*/) const
  {
    Assert(false, ExcNotImplemented());
    return Tensor<1,2>();
  }
  
  
  void
  LSingularityGradFunction::gradient_list (
    const std::vector<Point<2> >& /*points*/,
    std::vector<Tensor<1,2> >& /*gradients*/,
    const unsigned int /*component*/) const
  {
    Assert(false, ExcNotImplemented());
  }
  
  
  void
  LSingularityGradFunction::vector_gradient_list (
    const std::vector<Point<2> >& /*points*/,
    std::vector<std::vector<Tensor<1,2> > >& /*gradients*/) const
  {
    Assert(false, ExcNotImplemented());
  }
  
//////////////////////////////////////////////////////////////////////
  
  template <int dim>
  double
  SlitSingularityFunction<dim>::value (
    const Point<dim>   &p,
    const unsigned int) const
  {
    double x = p(0);
    double y = p(1);
    
    double phi = std::atan2(x,y)+M_PI;
    double r2 = x*x+y*y;
    
    return std::pow(r2,.25) * std::sin(.5*phi);
  }
  
  
  template <int dim>
  void
  SlitSingularityFunction<dim>::value_list (
    const std::vector<Point<dim> > &points,
    std::vector<double>            &values,
    const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	double x = points[i](0);
	double y = points[i](1);
	
	double phi = std::atan2(x,y)+M_PI;
	double r2 = x*x+y*y;
	
	values[i] = std::pow(r2,.25) * std::sin(.5*phi);
      }
  }
  
  
  template <int dim>
  void
  SlitSingularityFunction<dim>::vector_value_list (
    const std::vector<Point<dim> >& points,
    std::vector<Vector<double> >& values) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	Assert (values[i].size() == 1,
		ExcDimensionMismatch(values[i].size(), 1));
	
	double x = points[i](0);
	double y = points[i](1);
	
	double phi = std::atan2(x,y)+M_PI;
	double r2 = x*x+y*y;
	
	values[i](0) = std::pow(r2,.25) * std::sin(.5*phi);
      }
  }
  
  
  template <int dim>
  double
  SlitSingularityFunction<dim>::laplacian (const Point<dim>   &,
					   const unsigned int) const
  {
    return 0.;
  }
  
  
  template <int dim>
  void
  SlitSingularityFunction<dim>::laplacian_list (
    const std::vector<Point<dim> > &points,
    std::vector<double>            &values,
    const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      values[i] = 0.;
  }
  
  
  template <int dim>
  Tensor<1,dim>
  SlitSingularityFunction<dim>::gradient (const Point<dim>   &p,
					  const unsigned int) const
  {
    double x = p(0);
    double y = p(1);
    double phi = std::atan2(x,y)+M_PI;
    double r64 = std::pow(x*x+y*y,3./4.);
    
    Tensor<1,dim> result;
    result[0] = 1./2.*(std::sin(1./2.*phi)*x + std::cos(1./2.*phi)*y)/r64;
    result[1] = 1./2.*(std::sin(1./2.*phi)*y - std::cos(1./2.*phi)*x)/r64;
    return result;
  }
  
  
  template <int dim>
  void
  SlitSingularityFunction<dim>::gradient_list (const std::vector<Point<dim> > &points,
					       std::vector<Tensor<1,dim> >    &gradients,
					       const unsigned int) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	double x = p(0);
	double y = p(1);
	double phi = std::atan2(x,y)+M_PI;
	double r64 = std::pow(x*x+y*y,3./4.);
	
	gradients[i][0] = 1./2.*(std::sin(1./2.*phi)*x + std::cos(1./2.*phi)*y)/r64;
	gradients[i][1] = 1./2.*(std::sin(1./2.*phi)*y - std::cos(1./2.*phi)*x)/r64;
	for (unsigned int d=2;d<dim;++d)
	  gradients[i][d] = 0.;
      }
  }
  
  template <int dim>
  void
  SlitSingularityFunction<dim>::vector_gradient_list (
    const std::vector<Point<dim> >& points,
    std::vector<std::vector<Tensor<1,dim> > >& gradients) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	Assert(gradients[i].size() == 1,
	       ExcDimensionMismatch(gradients[i].size(), 1));
	
	const Point<dim>& p = points[i];
	double x = p(0);
	double y = p(1);
	double phi = std::atan2(x,y)+M_PI;
	double r64 = std::pow(x*x+y*y,3./4.);
	
	gradients[i][0][0] = 1./2.*(std::sin(1./2.*phi)*x + std::cos(1./2.*phi)*y)/r64;
	gradients[i][0][1] = 1./2.*(std::sin(1./2.*phi)*y - std::cos(1./2.*phi)*x)/r64;
	for (unsigned int d=2;d<dim;++d)
	  gradients[i][0][d] = 0.;
      }
  }
  
//////////////////////////////////////////////////////////////////////
  
  
  double
  SlitHyperSingularityFunction::value (const Point<2>   &p,
				       const unsigned int) const
  {
    double x = p(0);
    double y = p(1);
    
    double phi = std::atan2(x,y)+M_PI;
    double r2 = x*x+y*y;
    
    return std::pow(r2,.125) * std::sin(.25*phi);
  }
  
  
  void
  SlitHyperSingularityFunction::value_list (
    const std::vector<Point<2> > &points,
    std::vector<double>            &values,
    const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	double x = points[i](0);
	double y = points[i](1);
	
	double phi = std::atan2(x,y)+M_PI;
	double r2 = x*x+y*y;
	
	values[i] = std::pow(r2,.125) * std::sin(.25*phi);
      }
  }
  
  
  void
  SlitHyperSingularityFunction::vector_value_list (
    const std::vector<Point<2> >& points,
    std::vector<Vector<double> >& values) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	Assert(values[i].size() == 1,
	       ExcDimensionMismatch(values[i].size(), 1));
	
	double x = points[i](0);
	double y = points[i](1);
	
	double phi = std::atan2(x,y)+M_PI;
	double r2 = x*x+y*y;
	
	values[i](0) = std::pow(r2,.125) * std::sin(.25*phi);
      }
  }
  
  
  double
  SlitHyperSingularityFunction::laplacian (
    const Point<2>   &,
    const unsigned int) const
  {
    return 0.;
  }
  
  
  void
  SlitHyperSingularityFunction::laplacian_list (
    const std::vector<Point<2> > &points,
    std::vector<double>            &values,
    const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      values[i] = 0.;
  }
  
  
  Tensor<1,2>
  SlitHyperSingularityFunction::gradient (
    const Point<2>   &p,
    const unsigned int) const
  {
    double x = p(0);
    double y = p(1);
    double phi = std::atan2(x,y)+M_PI;
    double r78 = std::pow(x*x+y*y,7./8.);
    
    
    Tensor<1,2> result;
    result[0] = 1./4.*(std::sin(1./4.*phi)*x + std::cos(1./4.*phi)*y)/r78;
    result[1] = 1./4.*(std::sin(1./4.*phi)*y - std::cos(1./4.*phi)*x)/r78;
    return result;
  }
  
  
  void
  SlitHyperSingularityFunction::gradient_list (
    const std::vector<Point<2> > &points,
    std::vector<Tensor<1,2> >    &gradients,
    const unsigned int) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<2>& p = points[i];
	double x = p(0);
	double y = p(1);
	double phi = std::atan2(x,y)+M_PI;
	double r78 = std::pow(x*x+y*y,7./8.);
	
	gradients[i][0] = 1./4.*(std::sin(1./4.*phi)*x + std::cos(1./4.*phi)*y)/r78;
	gradients[i][1] = 1./4.*(std::sin(1./4.*phi)*y - std::cos(1./4.*phi)*x)/r78;
      }
  }
  
  
  void
  SlitHyperSingularityFunction::vector_gradient_list (
    const std::vector<Point<2> > &points,
    std::vector<std::vector<Tensor<1,2> > >& gradients) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	Assert(gradients[i].size() == 1,
	       ExcDimensionMismatch(gradients[i].size(), 1));
	
	const Point<2>& p = points[i];
	double x = p(0);
	double y = p(1);
	double phi = std::atan2(x,y)+M_PI;
	double r78 = std::pow(x*x+y*y,7./8.);
	
	gradients[i][0][0] = 1./4.*(std::sin(1./4.*phi)*x + std::cos(1./4.*phi)*y)/r78;
	gradients[i][0][1] = 1./4.*(std::sin(1./4.*phi)*y - std::cos(1./4.*phi)*x)/r78;
      }
  }
  
//////////////////////////////////////////////////////////////////////
  
  template<int dim>
  JumpFunction<dim>::JumpFunction(const Point<dim> &direction,
				  const double      steepness)
		  :
		  direction(direction),
		  steepness(steepness)
  {
    switch (dim)
      {
	case 1:
	      angle = 0;
	      break;
	case 2:
	      angle = std::atan2(direction(0),direction(1));
	      break;
	case 3:
	      Assert(false, ExcNotImplemented());
      }
    sine = std::sin(angle);
    cosine = std::cos(angle);
  }
  
  
  
  template<int dim>
  double
  JumpFunction<dim>::value (const Point<dim>   &p,
			    const unsigned int) const
  {
    double x = steepness*(-cosine*p(0)+sine*p(1));
    return -std::atan(x);
  }
  
  
  
  template<int dim>
  void
  JumpFunction<dim>::value_list (const std::vector<Point<dim> > &p,
				 std::vector<double>          &values,
				 const unsigned int) const
  {
    Assert (values.size() == p.size(),
	    ExcDimensionMismatch(values.size(), p.size()));
    
    for (unsigned int i=0;i<p.size();++i)
      {
	double x = steepness*(-cosine*p[i](0)+sine*p[i](1));
	values[i] = -std::atan(x);
      }
  }
  
  
  template<int dim>
  double
  JumpFunction<dim>::laplacian (const Point<dim>   &p,
				const unsigned int) const
  {
    double x = steepness*(-cosine*p(0)+sine*p(1));
    double r = 1+x*x;
    return 2*steepness*steepness*x/(r*r);
  }
  
  
  template<int dim>
  void
  JumpFunction<dim>::laplacian_list (const std::vector<Point<dim> > &p,
				     std::vector<double>          &values,
				     const unsigned int) const
  {
    Assert (values.size() == p.size(),
	    ExcDimensionMismatch(values.size(), p.size()));
    
    double f = 2*steepness*steepness;
    
    for (unsigned int i=0;i<p.size();++i)
      {
	double x = steepness*(-cosine*p[i](0)+sine*p[i](1));
	double r = 1+x*x;
	values[i] = f*x/(r*r);
      }
  }
  
  
  
  template<int dim>
  Tensor<1,dim>
  JumpFunction<dim>::gradient (const Point<dim>   &p,
			       const unsigned int) const
  {
    double x = steepness*(-cosine*p(0)+sine*p(1));
    double r = -steepness*(1+x*x);
    Tensor<1,dim> erg;
    erg[0] = cosine*r;
    erg[1] = sine*r;
    return erg;
  }
  
  
  
  template<int dim>
  void
  JumpFunction<dim>::gradient_list (const std::vector<Point<dim> > &p,
				    std::vector<Tensor<1,dim> >  &gradients,
				    const unsigned int) const
  {
    Assert (gradients.size() == p.size(),
	    ExcDimensionMismatch(gradients.size(), p.size()));
    
    for (unsigned int i=0; i<p.size(); ++i)
      {
	double x = steepness*(cosine*p[i](0)+sine*p[i](1));
	double r = -steepness*(1+x*x);
	gradients[i][0] = cosine*r;
	gradients[i][1] = sine*r;
      }
  }
  
  
  
  template <int dim>
  unsigned int
  JumpFunction<dim>::memory_consumption () const
  {
				     // only simple data elements, so
				     // use sizeof operator
    return sizeof (*this);
  }
  
  
  
  
  
/* ---------------------- FourierCosineFunction ----------------------- */
  
  
  template <int dim>
  FourierCosineFunction<dim>::
  FourierCosineFunction (const Point<dim> &fourier_coefficients)
		  :
		  Function<dim> (1),
		  fourier_coefficients (fourier_coefficients)
  {}
  
  
  
  template <int dim>
  double
  FourierCosineFunction<dim>::value (const Point<dim>   &p,
				     const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    return std::cos(fourier_coefficients * p);
  }
  
  
  
  template <int dim>
  Tensor<1,dim>
  FourierCosineFunction<dim>::gradient (const Point<dim>   &p,
					const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    return -fourier_coefficients * std::sin(fourier_coefficients * p);
  }
  
  
  
  template <int dim>
  double
  FourierCosineFunction<dim>::laplacian (const Point<dim>   &p,
					 const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    return fourier_coefficients.square() * (-std::cos(fourier_coefficients * p));
  }
  
  
  
  
/* ---------------------- FourierSineFunction ----------------------- */
  
  
  
  template <int dim>
  FourierSineFunction<dim>::
  FourierSineFunction (const Point<dim> &fourier_coefficients)
		  :
		  Function<dim> (1),
		  fourier_coefficients (fourier_coefficients)
  {}
  
  
  
  template <int dim>
  double
  FourierSineFunction<dim>::value (const Point<dim>   &p,
				   const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    return std::sin(fourier_coefficients * p);
  }
  
  
  
  template <int dim>
  Tensor<1,dim>
  FourierSineFunction<dim>::gradient (const Point<dim>   &p,
				      const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    return fourier_coefficients * std::cos(fourier_coefficients * p);
  }
  
  
  
  template <int dim>
  double
  FourierSineFunction<dim>::laplacian (const Point<dim>   &p,
				       const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    return fourier_coefficients.square() * (-std::sin(fourier_coefficients * p));
  }
  



/* ---------------------- FourierSineSum ----------------------- */
  
  
  
  template <int dim>
  FourierSineSum<dim>::
  FourierSineSum (const std::vector<Point<dim> > &fourier_coefficients,
		  const std::vector<double>      &weights)
		  :
		  Function<dim> (1),
		  fourier_coefficients (fourier_coefficients),
		  weights (weights)
  {
    Assert (fourier_coefficients.size() > 0, ExcZero());
    Assert (fourier_coefficients.size() == weights.size(),
	    ExcDimensionMismatch(fourier_coefficients.size(),
				 weights.size()));
  }
  
  
  
  template <int dim>
  double
  FourierSineSum<dim>::value (const Point<dim>   &p,
			      const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));

    const unsigned int n = weights.size();
    double sum = 0;
    for (unsigned int s=0; s<n; ++s)
      sum += weights[s] * std::sin(fourier_coefficients[s] * p);
    
    return sum;
  }
  
  
  
  template <int dim>
  Tensor<1,dim>
  FourierSineSum<dim>::gradient (const Point<dim>   &p,
				 const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));

    const unsigned int n = weights.size();
    Tensor<1,dim> sum;
    for (unsigned int s=0; s<n; ++s)
      sum += fourier_coefficients[s] * std::cos(fourier_coefficients[s] * p);
    
    return sum;
  }
  
  
  
  template <int dim>
  double
  FourierSineSum<dim>::laplacian (const Point<dim>   &p,
				  const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));

    const unsigned int n = weights.size();
    double sum = 0;
    for (unsigned int s=0; s<n; ++s)
      sum -= fourier_coefficients[s].square() * std::sin(fourier_coefficients[s] * p);
    
    return sum;
  }



/* ---------------------- FourierCosineSum ----------------------- */
  
  
  
  template <int dim>
  FourierCosineSum<dim>::
  FourierCosineSum (const std::vector<Point<dim> > &fourier_coefficients,
		    const std::vector<double>      &weights)
		  :
		  Function<dim> (1),
		  fourier_coefficients (fourier_coefficients),
		  weights (weights)
  {
    Assert (fourier_coefficients.size() > 0, ExcZero());
    Assert (fourier_coefficients.size() == weights.size(),
	    ExcDimensionMismatch(fourier_coefficients.size(),
				 weights.size()));
  }
  
  
  
  template <int dim>
  double
  FourierCosineSum<dim>::value (const Point<dim>   &p,
				const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));

    const unsigned int n = weights.size();
    double sum = 0;
    for (unsigned int s=0; s<n; ++s)
      sum += weights[s] * std::cos(fourier_coefficients[s] * p);
    
    return sum;
  }
  
  
  
  template <int dim>
  Tensor<1,dim>
  FourierCosineSum<dim>::gradient (const Point<dim>   &p,
				   const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));

    const unsigned int n = weights.size();
    Tensor<1,dim> sum;
    for (unsigned int s=0; s<n; ++s)
      sum -= fourier_coefficients[s] * std::sin(fourier_coefficients[s] * p);
    
    return sum;
  }
  
  
  
  template <int dim>
  double
  FourierCosineSum<dim>::laplacian (const Point<dim>   &p,
				    const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));

    const unsigned int n = weights.size();
    double sum = 0;
    for (unsigned int s=0; s<n; ++s)
      sum -= fourier_coefficients[s].square() * std::cos(fourier_coefficients[s] * p);
    
    return sum;
  }




/* ---------------------- Monomial ----------------------- */
  
  
  
  template <int dim>
  Monomial<dim>::
  Monomial (const Tensor<1,dim> &exponents,
	    const unsigned int   n_components)
		  :
		  Function<dim> (n_components),
		  exponents (exponents)
  {}
  
  
  
  template <int dim>
  double
  Monomial<dim>::value (const Point<dim>   &p,
			const unsigned int  component) const
  {
    Assert (component<this->n_components,
	    ExcIndexRange(component, 0, this->n_components)) ;

    double prod = 1;
    for (unsigned int s=0; s<dim; ++s)
      prod *= std::pow(p[s], exponents[s]);
    
    return prod;
  }



  template <int dim>
  void
  Monomial<dim>::vector_value (const Point<dim>   &p,
			       Vector<double>     &values) const
  {
    Assert (values.size() == this->n_components,
	    ExcDimensionMismatch (values.size(), this->n_components));

    for (unsigned int i=0; i<values.size(); ++i)
      values(i) = Monomial<dim>::value(p,i);
  }
  
  
  
  template <int dim>
  Tensor<1,dim>
  Monomial<dim>::gradient (const Point<dim>   &p,
			   const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1)) ;

    Tensor<1,dim> r;
    for (unsigned int d=0; d<dim; ++d)
      {
	double prod = 1;
	for (unsigned int s=0; s<dim; ++s)
	  prod *= (s==d
		   ?
		   exponents[s] * std::pow(p[s], exponents[s]-1)
		   :
		   std::pow(p[s], exponents[s]));

	r[d] = prod;
      }
  
    return r;
  }



  template<int dim>
  void
  Monomial<dim>::value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int i=0; i<points.size(); ++i)
      values[i] = Monomial<dim>::value (points[i], component);
  }
  
  
  template <int dim>
  Bessel1<dim>::Bessel1(
    const unsigned int order,
    const double wave_number,
    const Point<dim> center)
		  :
		  order(order),
		  wave_number(wave_number),
		  center(center)
  {}
  
  template <int dim>
  double
  Bessel1<dim>::value(const Point<dim>& p, const unsigned int) const
  {
    Assert(dim==2, ExcNotImplemented());
    const double r = p.distance(center);
#ifdef HAVE_JN_
    return jn(order, r*wave_number);
#else
    Assert(false, ExcMessage("Bessel function jn was not found by configure"));
    return r;
#endif
  }
  
  
  template <int dim>
  void
  Bessel1<dim>::value_list (
    const std::vector<Point<dim> > &points,
    std::vector<double>            &values,
    const unsigned int) const
  {
    Assert(dim==2, ExcNotImplemented());
    AssertDimension(points.size(), values.size());
    for (unsigned int k=0;k<points.size();++k)
      {
#ifdef HAVE_JN_
	const double r = points[k].distance(center);
	values[k] = jn(order, r*wave_number);
#else
    Assert(false, ExcMessage("Bessel function jn was not found by configure"));
#endif    
      }
  }

  
  template <int dim>
  Tensor<1,dim>
  Bessel1<dim>::gradient (const Point<dim>   &p,
			  const unsigned int) const
  {
    Assert(dim==2, ExcNotImplemented());
    const double r = p.distance(center);
    const double co = (r==0.) ? 0. : (p(0)-center(0))/r;
    const double si = (r==0.) ? 0. : (p(1)-center(1))/r;
    
    const double dJn = (order==0)
		       ? (-jn(1, r*wave_number))
		       : (.5*(jn(order-1, wave_number*r) -jn(order+1, wave_number*r)));
    Tensor<1,dim> result;
    result[0] = wave_number * co * dJn;
    result[1] = wave_number * si * dJn;
#ifdef HAVE_JN_
    return result;
#else
    Assert(false, ExcMessage("Bessel function jn was not found by configure"));
    return Tensor<1,dim>();
#endif
  }
  

  
  template <int dim>
  void
  Bessel1<dim>::gradient_list (
    const std::vector<Point<dim> > &points,
    std::vector<Tensor<1,dim> >    &gradients,
    const unsigned int) const
  {
    Assert(dim==2, ExcNotImplemented());
    AssertDimension(points.size(), gradients.size());
    for (unsigned int k=0;k<points.size();++k)
      {
	const Point<dim>& p = points[k];
	const double r = p.distance(center);
	const double co = (r==0.) ? 0. : (p(0)-center(0))/r;
	const double si = (r==0.) ? 0. : (p(1)-center(1))/r;
	
#ifdef HAVE_JN_
	const double dJn = (order==0)
			   ? (-jn(1, r*wave_number))
			   : (.5*(jn(order-1, wave_number*r) -jn(order+1, wave_number*r)));
#else
	const double dJn = 0.;
	Assert(false, ExcMessage("Bessel function jn was not found by configure"));
#endif
	Tensor<1,dim>& result = gradients[k];
	result[0] = wave_number * co * dJn;
	result[1] = wave_number * si * dJn;
      }
  }
  

// explicit instantiations  
  template class SquareFunction<1>;
  template class SquareFunction<2>;
  template class SquareFunction<3>;
  template class Q1WedgeFunction<1>;
  template class Q1WedgeFunction<2>;
  template class Q1WedgeFunction<3>;
  template class PillowFunction<1>;
  template class PillowFunction<2>;
  template class PillowFunction<3>;
  template class CosineFunction<1>;
  template class CosineFunction<2>;
  template class CosineFunction<3>;
  template class CosineGradFunction<1>;
  template class CosineGradFunction<2>;
  template class CosineGradFunction<3>;
  template class ExpFunction<1>;
  template class ExpFunction<2>;
  template class ExpFunction<3>;
  template class JumpFunction<1>;
  template class JumpFunction<2>;
  template class JumpFunction<3>;
  template class FourierCosineFunction<1>;
  template class FourierCosineFunction<2>;
  template class FourierCosineFunction<3>;
  template class FourierSineFunction<1>;
  template class FourierSineFunction<2>;
  template class FourierSineFunction<3>;
  template class FourierCosineSum<1>;
  template class FourierCosineSum<2>;
  template class FourierCosineSum<3>;
  template class FourierSineSum<1>;
  template class FourierSineSum<2>;
  template class FourierSineSum<3>;
  template class SlitSingularityFunction<2>;
  template class SlitSingularityFunction<3>;
  template class Monomial<1>;
  template class Monomial<2>;
  template class Monomial<3>;
  template class Bessel1<2>;
  template class Bessel1<3>;
}

DEAL_II_NAMESPACE_CLOSE
