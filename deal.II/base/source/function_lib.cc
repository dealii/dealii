// $Id$
// (c) Guido Kanschat, 1999

#include <base/point.h>
#include <base/function_lib.h>

//TODO: Derivatives in 3d wrong

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

template PillowFunction<1>;
template PillowFunction<2>;
template PillowFunction<3>;

