//----------------------------  function_lib.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function_lib.cc  ---------------------------


#include <base/tensor.h>
#include <base/point.h>
#include <base/function_lib.h>

#include <cmath>


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
  
  
// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
  using namespace std;
#endif
  
  template<int dim>
  void
  SquareFunction<dim>::value_list (const typename std::vector<Point<dim> > &points,
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
  SquareFunction<dim>::laplacian_list (const typename std::vector<Point<dim> > &points,
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
    return p*2;
  }
  
  
  
  template<int dim>
  void
  SquareFunction<dim>::gradient_list (const typename std::vector<Point<dim> > &points,
				      typename std::vector<Tensor<1,dim> >    &gradients,
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
    return p(0)*p(1);
  }
  
  
  
  template<int dim>
  void
  Q1WedgeFunction<dim>::value_list (const typename std::vector<Point<dim> > &points,
				    std::vector<double>            &values,
				    const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      {
	const Point<dim>& p = points[i];
	values[i] = p(0)*p(1);
      }
  }
  
  
  template<int dim>
  double
  Q1WedgeFunction<dim>::laplacian (const Point<dim>   &,
				   const unsigned int) const
  {
    return 0.;
  }
  
  
  template<int dim>
  void
  Q1WedgeFunction<dim>::laplacian_list (const typename std::vector<Point<dim> > &points,
					std::vector<double>            &values,
					const unsigned int) const
  {
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
    Tensor<1,dim> erg;
    erg[0] = p(1);
    erg[1] = p(0);
    return erg;
  }
  
  
  
  template<int dim>
  void
  Q1WedgeFunction<dim>::gradient_list (const typename std::vector<Point<dim> > &points,
				       typename std::vector<Tensor<1,dim> >    &gradients,
				       const unsigned int) const
  {
    Assert (gradients.size() == points.size(),
	    ExcDimensionMismatch(gradients.size(), points.size()));
    
    for (unsigned int i=0; i<points.size(); ++i)
      {
	gradients[i][0] = points[i](1);
	gradients[i][1] = points[i](0);
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
  PillowFunction<dim>::value_list (const typename std::vector<Point<dim> > &points,
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
  PillowFunction<dim>::laplacian_list (const typename std::vector<Point<dim> > &points,
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
  PillowFunction<dim>::gradient_list (const typename std::vector<Point<dim> > &points,
				      typename std::vector<Tensor<1,dim> >    &gradients,
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
  CosineFunction<dim>::value_list (const typename std::vector<Point<dim> > &points,
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
		  values[i] = std::cos(M_PI_2*p(0));
		  break;
	    case 2:
		  values[i] = std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
		  break;
	    case 3:
		  values[i] = std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
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
  CosineFunction<dim>::laplacian_list (const typename std::vector<Point<dim> > &points,
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
		  values[i] = -M_PI_2*M_PI_2* std::cos(M_PI_2*p(0));
		  break;
	    case 2:
		  values[i] = -2*M_PI_2*M_PI_2* std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1));
		  break;
	    case 3:
		  values[i] = -3*M_PI_2*M_PI_2* std::cos(M_PI_2*p(0)) * std::cos(M_PI_2*p(1)) * std::cos(M_PI_2*p(2));
		  break;
	    default:
		  Assert(false, ExcNotImplemented());
	  }
      }
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
  CosineFunction<dim>::gradient_list (const typename std::vector<Point<dim> > &points,
				      typename std::vector<Tensor<1,dim> >    &gradients,
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
  CosineFunction<dim>::hessian_list (const typename std::vector<Point<dim> > &points,
				     typename std::vector<Tensor<2,dim> >    &hessians,
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
  ExpFunction<dim>::value_list (const typename std::vector<Point<dim> > &points,
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
  ExpFunction<dim>::laplacian_list (const typename std::vector<Point<dim> > &points,
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
  ExpFunction<dim>::gradient_list (const typename std::vector<Point<dim> > &points,
				   typename std::vector<Tensor<1,dim> >    &gradients,
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
  
//////////////////////////////////////////////////////////////////////
  
  
  double
  SlitSingularityFunction::value (const Point<2>   &p,
				  const unsigned int) const
  {
    double x = p(0);
    double y = p(1);
    
    double phi = std::atan2(x,y)+M_PI;
    double r2 = x*x+y*y;
    
    return std::pow(r2,.25) * std::sin(.5*phi);
  }
  
  
  void
  SlitSingularityFunction::value_list (const std::vector<Point<2> > &points,
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
  
  
  double
  SlitSingularityFunction::laplacian (const Point<2>   &,
				      const unsigned int) const
  {
    return 0.;
  }
  
  
  void
  SlitSingularityFunction::laplacian_list (const std::vector<Point<2> > &points,
					   std::vector<double>            &values,
					   const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<points.size();++i)
      values[i] = 0.;
  }
  
  
  Tensor<1,2>
  SlitSingularityFunction::gradient (const Point<2>   &p,
				     const unsigned int) const
  {
    double x = p(0);
    double y = p(1);
    double phi = std::atan2(x,y)+M_PI;
    double r64 = std::pow(x*x+y*y,3./4.);
    
    
    Tensor<1,2> result;
    result[0] = 1./2.*(std::sin(1./2.*phi)*x + std::cos(1./2.*phi)*y)/r64;
    result[1] = 1./2.*(std::sin(1./2.*phi)*y - std::cos(1./2.*phi)*x)/r64;
    return result;
  }
  
  
  void
  SlitSingularityFunction::gradient_list (const std::vector<Point<2> > &points,
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
	double r64 = std::pow(x*x+y*y,3./4.);
	
	gradients[i][0] = 1./2.*(std::sin(1./2.*phi)*x + std::cos(1./2.*phi)*y)/r64;
	gradients[i][1] = 1./2.*(std::sin(1./2.*phi)*y - std::cos(1./2.*phi)*x)/r64;
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
  JumpFunction<dim>::value_list (const typename std::vector<Point<dim> > &p,
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
  JumpFunction<dim>::laplacian_list (const typename std::vector<Point<dim> > &p,
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
  JumpFunction<dim>::gradient_list (const typename std::vector<Point<dim> > &p,
				    typename std::vector<Tensor<1,dim> >  &gradients,
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
  };
  
  
  
  
  
/* ---------------------- FourierSineFunction ----------------------- */
  
  
  template <int dim>
  FourierCosineFunction<dim>::
  FourierCosineFunction (const Point<dim> &fourier_coefficients)
		  :
		  Function<dim> (1),
    fourier_coefficients (fourier_coefficients)
  {};
  
  
  
  template <int dim>
  double
  FourierCosineFunction<dim>::value (const Point<dim>   &p,
				     const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    double val=1;
    for (unsigned int i=0; i<dim; ++i)
      val *= std::cos(fourier_coefficients[i] * p[i]);
    return val;
  };
  
  
  
  template <int dim>
  Tensor<1,dim>
  FourierCosineFunction<dim>::gradient (const Point<dim>   &p,
					const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    Tensor<1,dim> grad;
    for (unsigned int i=0; i<dim; ++i)
      grad[i] = 1;
    
    for (unsigned int i=0; i<dim; ++i)
      {
	const double cos_i = std::cos(fourier_coefficients[i] * p[i]);
	const double sin_i = std::sin(fourier_coefficients[i] * p[i]);
	
	for (unsigned int d=0; d<dim; ++d)
	  if (d==i)
	    grad[d] *= - fourier_coefficients[i] * sin_i;
	  else
	    grad[d] *= cos_i;
      };
    
    return grad;
  };
  
  
  
  template <int dim>
  double
  FourierCosineFunction<dim>::laplacian (const Point<dim>   &p,
					 const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    double val = -(fourier_coefficients*fourier_coefficients);
    for (unsigned int i=0; i<dim; ++i)
      val *= std::cos(fourier_coefficients[i] * p[i]);
    return val;
  };
  
  
  
  
/* ---------------------- FourierSineFunction ----------------------- */
  
  
  
  template <int dim>
  FourierSineFunction<dim>::
  FourierSineFunction (const Point<dim> &fourier_coefficients)
		  :
		  Function<dim> (1),
    fourier_coefficients (fourier_coefficients)
  {};
  
  
  
  template <int dim>
  double
  FourierSineFunction<dim>::value (const Point<dim>   &p,
				   const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    double val=1;
    for (unsigned int i=0; i<dim; ++i)
      val *= std::sin(fourier_coefficients[i] * p[i]);
    return val;
  };
  
  
  
  template <int dim>
  Tensor<1,dim>
  FourierSineFunction<dim>::gradient (const Point<dim>   &p,
				      const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    Tensor<1,dim> grad;
    for (unsigned int i=0; i<dim; ++i)
      grad[i] = 1;
    
    for (unsigned int i=0; i<dim; ++i)
      {
	const double cos_i = std::cos(fourier_coefficients[i] * p[i]);
	const double sin_i = std::sin(fourier_coefficients[i] * p[i]);
	
	for (unsigned int d=0; d<dim; ++d)
	  if (d==i)
	    grad[d] *= fourier_coefficients[i] * cos_i;
	  else
	    grad[d] *= sin_i;
      };
    
    return grad;
  };
  
  
  
  template <int dim>
  double
  FourierSineFunction<dim>::laplacian (const Point<dim>   &p,
				       const unsigned int  component) const
  {
    Assert (component==0, ExcIndexRange(component,0,1));
    double val = -(fourier_coefficients*fourier_coefficients);
    for (unsigned int i=0; i<dim; ++i)
      val *= std::sin(fourier_coefficients[i] * p[i]);
    return val;
  };
  
  
  
  
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
  
  
};
