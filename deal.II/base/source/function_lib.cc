//----------------------------  function_lib.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function_lib.cc  ---------------------------


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



template<int dim>
double
SquareFunction<dim>::value (const Point<dim>   &p,
			    const unsigned int) const
{
  return p.square();
}



template<int dim>
void
SquareFunction<dim>::value_list (const vector<Point<dim> > &points,
				 vector<double>            &values,
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
SquareFunction<dim>::laplacian_list (const vector<Point<dim> > &points,
				     vector<double>            &values,
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
SquareFunction<dim>::gradient_list (const vector<Point<dim> > &points,
				    vector<Tensor<1,dim> >    &gradients,
				    const unsigned int) const
{
  Assert (gradients.size() == points.size(),
	  ExcDimensionMismatch(gradients.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    gradients[i] = points[i]*2;
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
PillowFunction<dim>::value_list (const vector<Point<dim> > &points,
				 vector<double>            &values,
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
PillowFunction<dim>::laplacian_list (const vector<Point<dim> > &points,
				     vector<double>            &values,
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
PillowFunction<dim>::gradient_list (const vector<Point<dim> > &points,
				    vector<Tensor<1,dim> >    &gradients,
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
	    return cos(M_PI_2*p(0));
      case 2:
	    return cos(M_PI_2*p(0)) * cos(M_PI_2*p(1));
      case 3:
	    return cos(M_PI_2*p(0)) * cos(M_PI_2*p(1)) * cos(M_PI_2*p(2));
      default:
	    Assert(false, ExcNotImplemented());
    }
  return 0.;
}

template<int dim>
void
CosineFunction<dim>::value_list (const vector<Point<dim> > &points,
				 vector<double>            &values,
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
		values[i] = cos(M_PI_2*p(0));
		break;
	  case 2:
		values[i] = cos(M_PI_2*p(0)) * cos(M_PI_2*p(1));
		break;
	  case 3:
		values[i] = cos(M_PI_2*p(0)) * cos(M_PI_2*p(1)) * cos(M_PI_2*p(2));
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
	    return -M_PI_2*M_PI_2* cos(M_PI_2*p(0));
      case 2:
	    return -2*M_PI_2*M_PI_2* cos(M_PI_2*p(0)) * cos(M_PI_2*p(1));
      case 3:
	    return -3*M_PI_2*M_PI_2* cos(M_PI_2*p(0)) * cos(M_PI_2*p(1)) * cos(M_PI_2*p(2));
      default:
	    Assert(false, ExcNotImplemented());
    }
  return 0.;
}

template<int dim>
void
CosineFunction<dim>::laplacian_list (const vector<Point<dim> > &points,
				     vector<double>            &values,
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
		values[i] = -M_PI_2*M_PI_2* cos(M_PI_2*p(0));
		break;
	  case 2:
		values[i] = -2*M_PI_2*M_PI_2* cos(M_PI_2*p(0)) * cos(M_PI_2*p(1));
		break;
	  case 3:
		values[i] = -3*M_PI_2*M_PI_2* cos(M_PI_2*p(0)) * cos(M_PI_2*p(1)) * cos(M_PI_2*p(2));
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
	    result[0] = -M_PI_2* sin(M_PI_2*p(0));
	    break;
      case 2:
	    result[0] = -M_PI_2* sin(M_PI_2*p(0)) * cos(M_PI_2*p(1));
	    result[1] = -M_PI_2* cos(M_PI_2*p(0)) * sin(M_PI_2*p(1));
	    break;
      case 3:
	    result[0] = -M_PI_2* sin(M_PI_2*p(0)) * cos(M_PI_2*p(1)) * cos(M_PI_2*p(2));
	    result[1] = -M_PI_2* cos(M_PI_2*p(0)) * sin(M_PI_2*p(1)) * cos(M_PI_2*p(2));
	    result[2] = -M_PI_2* cos(M_PI_2*p(0)) * cos(M_PI_2*p(1)) * sin(M_PI_2*p(2));
	    break;
      default:
	    Assert(false, ExcNotImplemented());
    }
  return result;
}

template<int dim>
void
CosineFunction<dim>::gradient_list (const vector<Point<dim> > &points,
				    vector<Tensor<1,dim> >    &gradients,
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
		gradients[i][0] = -M_PI_2* sin(M_PI_2*p(0));
		break;
	  case 2:
		gradients[i][0] = -M_PI_2* sin(M_PI_2*p(0)) * cos(M_PI_2*p(1));
		gradients[i][1] = -M_PI_2* cos(M_PI_2*p(0)) * sin(M_PI_2*p(1));
		break;
	  case 3:
		gradients[i][0] = -M_PI_2* sin(M_PI_2*p(0)) * cos(M_PI_2*p(1)) * cos(M_PI_2*p(2));
		gradients[i][1] = -M_PI_2* cos(M_PI_2*p(0)) * sin(M_PI_2*p(1)) * cos(M_PI_2*p(2));
		gradients[i][2] = -M_PI_2* cos(M_PI_2*p(0)) * cos(M_PI_2*p(1)) * sin(M_PI_2*p(2));
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
	    return exp(p(0));
      case 2:
	    return exp(p(0)) * exp(p(1));
      case 3:
	    return exp(p(0)) * exp(p(1)) * exp(p(2));
      default:
	    Assert(false, ExcNotImplemented());
    }
  return 0.;
}

template<int dim>
void
ExpFunction<dim>::value_list (const vector<Point<dim> > &points,
			      vector<double>            &values,
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
		values[i] = exp(p(0));
		break;
	  case 2:
		values[i] = exp(p(0)) * exp(p(1));
		break;
	  case 3:
		values[i] = exp(p(0)) * exp(p(1)) * exp(p(2));
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
	    return exp(p(0));
      case 2:
	    return 2 * exp(p(0)) * exp(p(1));
      case 3:
	    return 3 * exp(p(0)) * exp(p(1)) * exp(p(2));
      default:
	    Assert(false, ExcNotImplemented());
    }
  return 0.;
}

template<int dim>
void
ExpFunction<dim>::laplacian_list (const vector<Point<dim> > &points,
				  vector<double>            &values,
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
		values[i] = exp(p(0));
		break;
	  case 2:
		values[i] = 2 * exp(p(0)) * exp(p(1));
		break;
	  case 3:
		values[i] = 3 * exp(p(0)) * exp(p(1)) * exp(p(2));
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
	    result[0] = exp(p(0));
	    break;
      case 2:
	    result[0] = exp(p(0)) * exp(p(1));
	    result[1] = exp(p(0)) * exp(p(1));
	    break;
      case 3:
	    result[0] = exp(p(0)) * exp(p(1)) * exp(p(2));
	    result[1] = exp(p(0)) * exp(p(1)) * exp(p(2));
	    result[2] = exp(p(0)) * exp(p(1)) * exp(p(2));
	    break;
      default:
	    Assert(false, ExcNotImplemented());
    }
  return result;
}

template<int dim>
void
ExpFunction<dim>::gradient_list (const vector<Point<dim> > &points,
				 vector<Tensor<1,dim> >    &gradients,
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
		gradients[i][0] = exp(p(0));
		break;
	  case 2:
		gradients[i][0] = exp(p(0)) * exp(p(1));
		gradients[i][1] = exp(p(0)) * exp(p(1));
		break;
	  case 3:
		gradients[i][0] = exp(p(0)) * exp(p(1)) * exp(p(2));
		gradients[i][1] = exp(p(0)) * exp(p(1)) * exp(p(2));
		gradients[i][2] = exp(p(0)) * exp(p(1)) * exp(p(2));
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
  
  double phi = atan2(y,-x)+M_PI;
  double r2 = x*x+y*y;

  return pow(r2,1./3.) * sin(2./3.*phi);
}


void
LSingularityFunction::value_list (const vector<Point<2> > &points,
				  vector<double>            &values,
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
	  double phi = atan2(y,-x)+M_PI;
	  double r2 = x*x+y*y;

	  values[i] = pow(r2,1./3.) * sin(2./3.*phi);
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
LSingularityFunction::laplacian_list (const vector<Point<2> > &points,
				      vector<double>            &values,
				      const unsigned int) const
{
  Assert (values.size() == points.size(),
	  ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0;i<points.size();++i)
    values[i] = 0.;
}


//TODO: Implement derivatives

Tensor<1,2>
LSingularityFunction::gradient (const Point<2>   &p,
				const unsigned int) const
{
  Assert(false, ExcNotImplemented());
  double x = p(0);
  double y = p(1);

  if ((x>=0) && (y>=0))
    {
//TODO: should return infinity, but how to do that?
      static const double infty[2] = {0., 0.};
      return Tensor<1,2>(infty);
    }
  
//  double phi = atan2(y,-x)+M_PI;
//  double r2 = x*x+y*y;

  return Tensor<1,2>();
}


void
LSingularityFunction::gradient_list (const vector<Point<2> > &points,
				     vector<Tensor<1,2> >    &gradients,
				     const unsigned int) const
{
  Assert (gradients.size() == points.size(),
	  ExcDimensionMismatch(gradients.size(), points.size()));
  Assert(false, ExcNotImplemented());
}

//////////////////////////////////////////////////////////////////////


double
SlitSingularityFunction::value (const Point<2>   &p,
				const unsigned int) const
{
  double x = p(0);
  double y = p(1);

  double phi = atan2(x,-y)+M_PI;
  double r2 = x*x+y*y;

  return pow(r2,.25) * sin(.5*phi);
}


void
SlitSingularityFunction::value_list (const vector<Point<2> > &points,
				     vector<double>            &values,
				     const unsigned int) const
{
  Assert (values.size() == points.size(),
	  ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0;i<points.size();++i)
    {
      double x = points[i](0);
      double y = points[i](1);

      double phi = atan2(x,-y)+M_PI;
      double r2 = x*x+y*y;

      values[i] = pow(r2,.25) * sin(.5*phi);
    }
}


double
SlitSingularityFunction::laplacian (const Point<2>   &,
				    const unsigned int) const
{
  return 0.;
}


void
SlitSingularityFunction::laplacian_list (const vector<Point<2> > &points,
					 vector<double>            &values,
					 const unsigned int) const
{
  Assert (values.size() == points.size(),
	  ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0;i<points.size();++i)
    values[i] = 0.;
}


//TODO: Implement derivatives

Tensor<1,2>
SlitSingularityFunction::gradient (const Point<2>   &/*p*/,
				   const unsigned int) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1,2>();
}


void
SlitSingularityFunction::gradient_list (const vector<Point<2> > &points,
					vector<Tensor<1,2> >    &gradients,
					const unsigned int) const
{
  Assert (gradients.size() == points.size(),
	  ExcDimensionMismatch(gradients.size(), points.size()));
  Assert(false, ExcNotImplemented());
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
	    angle = atan2(direction(0),direction(1));
	    break;
      case 3:
	    Assert(false, ExcNotImplemented());
    }
  sine = sin(angle);
  cosine = cos(angle);
}

	    

template<int dim>
double
JumpFunction<dim>::value (const Point<dim>   &p,
			  const unsigned int) const
{
  double x = steepness*(-cosine*p(0)+sine*p(1));
  return -atan(x);
}



template<int dim>
void
JumpFunction<dim>::value_list (const vector<Point<dim> > &p,
			       vector<double>          &values,
			       const unsigned int) const
{
  Assert (values.size() == p.size(),
	  ExcDimensionMismatch(values.size(), p.size()));

  for (unsigned int i=0;i<p.size();++i)
    {
      double x = steepness*(-cosine*p[i](0)+sine*p[i](1));
      values[i] = -atan(x);
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
JumpFunction<dim>::laplacian_list (const vector<Point<dim> > &p,
				   vector<double>          &values,
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
JumpFunction<dim>::gradient_list (const vector<Point<dim> > &p,
				  vector<Tensor<1,dim> >  &gradients,
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


template SquareFunction<1>;
template SquareFunction<2>;
template SquareFunction<3>;
template PillowFunction<1>;
template PillowFunction<2>;
template PillowFunction<3>;
template CosineFunction<1>;
template CosineFunction<2>;
template CosineFunction<3>;
template ExpFunction<1>;
template ExpFunction<2>;
template ExpFunction<3>;
template JumpFunction<1>;
template JumpFunction<2>;
template JumpFunction<3>;
