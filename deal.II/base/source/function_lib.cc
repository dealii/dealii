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

//TODO: Derivatives in 3d wrong (GK!)

template<int dim>
double
PillowFunction<dim>::value (const Point<dim>   &p,
			    const unsigned int) const
{
  switch(dim)
    {
      case 1:
	    return 1.-p(0)*p(0);
      case 2:
	    return (1.-p(0)*p(0))*(1.-p(1)*p(1));
      case 3:
	    return (1.-p(0)*p(0))*(1.-p(1)*p(1))*(1.-p(2)*p(2));
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
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0;i<points.size();++i)
    {
      const Point<dim>& p = points[i];
      switch(dim)
	{
	  case 1:
		values[i] = 1.-p(0)*p(0);
		break;
	  case 2:
		values[i] = (1.-p(0)*p(0))*(1.-p(1)*p(1));
		break;
	  case 3:
		values[i] = (1.-p(0)*p(0))*(1.-p(1)*p(1))*(1.-p(2)*p(2));
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
	    return 2.;
      case 2:
	    return 2.*((1.-p(0)*p(0))+(1.-p(1)*p(1)));
      case 3:
	    return 2.*((1.-p(0)*p(0))*(1.-p(1)*p(1))
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
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0;i<points.size();++i)
    {
      const Point<dim>& p = points[i];
      switch(dim)
	{
	  case 1:
		values[i] = 2.;
		break;
	  case 2:
		values[i] = 2.*((1.-p(0)*p(0))+(1.-p(1)*p(1)));
		break;
	  case 3:
		values[i] = 2.*((1.-p(0)*p(0))+(1.-p(1)*p(1))+(1.-p(2)*p(2)));
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
	    result[0] = 2.*p(0);
	    break;
      case 2:
	    result[0] = 2.*p(0)*(1.-p(1)*p(1));
	    result[1] = 2.*p(1)*(1.-p(0)*p(0));
	    break;
      case 3:
	    result[0] = 2.*p(0)*(1.-p(1)*p(1))*(1.-p(2)*p(2));
	    result[1] = 2.*p(1)*(1.-p(0)*p(0))*(1.-p(2)*p(2));
	    result[2] = 2.*p(2)*(1.-p(0)*p(0))*(1.-p(1)*p(1));
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
	  ExcVectorHasWrongSize(gradients.size(), points.size()));

  for (unsigned int i=0;i<points.size();++i)
    {
      const Point<dim>& p = points[i];
      switch(dim)
	{
	  case 1:
		gradients[i][0] = 2.*p(0);
		break;
	  case 2:
		gradients[i][0] = 2.*p(0)*(1.-p(1)*p(1));
		gradients[i][1] = 2.*p(1)*(1.-p(0)*p(0));
		return;
	  case 3:
		gradients[i][0] = 2.*p(0)*(1.-p(1)*p(1))*(1.-p(2)*p(2));
		gradients[i][1] = 2.*p(1)*(1.-p(0)*p(0))*(1.-p(2)*p(2));
		gradients[i][2] = 2.*p(2)*(1.-p(0)*p(0))*(1.-p(1)*p(1));
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
	  ExcVectorHasWrongSize(values.size(), points.size()));

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
	  ExcVectorHasWrongSize(values.size(), points.size()));

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
	  ExcVectorHasWrongSize(gradients.size(), points.size()));

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
		return;
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
	  ExcVectorHasWrongSize(values.size(), points.size()));

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
	  ExcVectorHasWrongSize(values.size(), points.size()));

  for (unsigned int i=0;i<points.size();++i)
    values[i] = 0.;
}


//TODO: Implement derivatives

Tensor<1,2>
LSingularityFunction::gradient (const Point<2>   &/*p*/,
			       const unsigned int) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1,2>();
}


void
LSingularityFunction::gradient_list (const vector<Point<2> > &points,
				    vector<Tensor<1,2> >    &gradients,
				    const unsigned int) const
{
  Assert (gradients.size() == points.size(),
	  ExcVectorHasWrongSize(gradients.size(), points.size()));
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
	  ExcVectorHasWrongSize(values.size(), points.size()));

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
	  ExcVectorHasWrongSize(values.size(), points.size()));

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
	  ExcVectorHasWrongSize(gradients.size(), points.size()));
  Assert(false, ExcNotImplemented());
}

//////////////////////////////////////////////////////////////////////

template PillowFunction<1>;
template PillowFunction<2>;
template PillowFunction<3>;
template CosineFunction<1>;
template CosineFunction<2>;
template CosineFunction<3>;
