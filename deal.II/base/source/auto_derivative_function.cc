//--------------------------------------------------------------------------
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
//--------------------------------------------------------------------------


#include <base/point.h>
#include <base/auto_derivative_function.h>
#include <lac/vector.h>

#include <cmath>


// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif


template <int dim>
AutoDerivativeFunction<dim>::AutoDerivativeFunction (const double hh,
						     const unsigned int n_components,
						     const double       initial_time):
		Function<dim>(n_components, initial_time),
		h(1),
		ht(dim),
		formula(Euler)
{
  set_h(hh);
  set_formula();
}


template <int dim>
AutoDerivativeFunction<dim>::~AutoDerivativeFunction ()
{}



template <int dim>
void
AutoDerivativeFunction<dim>::set_formula (const DifferenceFormula form)
{
  formula = form;
}


template <int dim>
void
AutoDerivativeFunction<dim>::set_h (const double hh)
{
  h=hh;
  for (unsigned int i=0; i<dim; ++i)
    ht[i][i]=h;
}


template <int dim>
Tensor<1,dim>
AutoDerivativeFunction<dim>::gradient (const Point<dim>   &p,
				       const unsigned int  comp) const
{
  Tensor<1,dim> grad;
  switch (formula)
    {
      case UpwindEuler:
      {
	Tensor<1,dim> q1;
	for (unsigned int i=0; i<dim; ++i)
	  {
	    q1=p-ht[i];
	    grad[i]=(value(p, comp)-value(q1, comp))/h;
	  }
	break;
      }
      case Euler:
      {
	Tensor<1,dim> q1, q2;
	for (unsigned int i=0; i<dim; ++i)
	  {
	    q1=p+ht[i];
	    q2=p-ht[i];
	    grad[i]=(value(q1, comp)-value(q2, comp))/(2*h);
	  }
	break;
      }
      case FourthOrder:
      {
	Tensor<1,dim> q1, q2, q3, q4;
	for (unsigned int i=0; i<dim; ++i)
	  {
	    q2=p+ht[i];
	    q1=q2+ht[i];
	    q3=p-ht[i];
	    q4=q3-ht[i];
	    grad[i]=(-  value(q1, comp)
		     +8*value(q2, comp)
		     -8*value(q3, comp)
		     +  value(q4, comp))/(12*h);
	    
	  }
	break;
      }
      default:
	    Assert(false, ExcInvalidFormula());
    }
  return grad;
}


template <int dim>
void AutoDerivativeFunction<dim>::vector_gradient (const Point<dim>       &p,
						   typename std::vector<Tensor<1,dim> > &gradients) const
{
  Assert (gradients.size() == n_components, ExcDimensionMismatch(gradients.size(), n_components));
  
  switch (formula)
    {
      case UpwindEuler:
      {
	Tensor<1,dim> q1;
	Vector<double> v(n_components), v1(n_components);
	const double h_inv=1./h;
	for (unsigned int i=0; i<dim; ++i)
	  {
	    q1=p-ht[i];
	    vector_value(p, v);
	    vector_value(q1, v1);
	    
	    for (unsigned int comp=0; comp<n_components; ++comp)
	      gradients[comp][i]=(v(comp)-v1(comp))*h_inv;
	  }
	break;
      }
      case Euler:
      {
	Tensor<1,dim> q1, q2;
	Vector<double> v1(n_components), v2(n_components);
	const double h_inv_2=1./(2*h);
	for (unsigned int i=0; i<dim; ++i)
	  {
	    q1=p+ht[i];
	    q2=p-ht[i];
	    vector_value(q1, v1);
	    vector_value(q2, v2);
	    
	    for (unsigned int comp=0; comp<n_components; ++comp)
	      gradients[comp][i]=(v1(comp)-v2(comp))*h_inv_2;
	  }
	break;
      }
      case FourthOrder:
      {
	Tensor<1,dim> q1, q2, q3, q4;
	Vector<double> v1(n_components), v2(n_components), v3(n_components), v4(n_components);
	const double h_inv_12=1./(12*h);
	for (unsigned int i=0; i<dim; ++i)
	  {
	    q2=p+ht[i];
	    q1=q2+ht[i];
	    q3=p-ht[i];
	    q4=q3-ht[i];
	    vector_value(q1, v1);
	    vector_value(q2, v2);
	    vector_value(q3, v3);
	    vector_value(q4, v4);
	    
	    for (unsigned int comp=0; comp<n_components; ++comp)
	      gradients[comp][i]=(-v1(comp)+8*v2(comp)-8*v3(comp)+v4(comp))*h_inv_12;
	  }
	break;
      }
      default:
	    Assert(false, ExcInvalidFormula());
    }
}


template <int dim>
void AutoDerivativeFunction<dim>::gradient_list (const typename std::vector<Point<dim> > &points,
						 typename std::vector<Tensor<1,dim> >    &gradients,
						 const unsigned int              comp) const
{
  Assert (gradients.size() == points.size(),
	  ExcDimensionMismatch(gradients.size(), points.size()));

  switch (formula)
    {
      case UpwindEuler:
      {
	Tensor<1,dim> q1;
	for (unsigned p=0; p<points.size(); ++p)
	  for (unsigned int i=0; i<dim; ++i)
	    {
	      q1=points[p]-ht[i];
	      gradients[p][i]=(value(points[p], comp)-value(q1, comp))/h;
	    }
	break;
      }
      case Euler:
      {
	Tensor<1,dim> q1, q2;
	for (unsigned p=0; p<points.size(); ++p)
	  for (unsigned int i=0; i<dim; ++i)
	    {
	      q1=points[p]+ht[i];
	      q2=points[p]-ht[i];
	      gradients[p][i]=(value(q1, comp)-value(q2, comp))/(2*h);
	    }
	break;
      }
      case FourthOrder:
      {
	Tensor<1,dim> q1, q2, q3, q4;
	for (unsigned p=0; p<points.size(); ++p)
	  for (unsigned int i=0; i<dim; ++i)
	    {
	      q2=points[p]+ht[i];
	      q1=q2+ht[i];
	      q3=points[p]-ht[i];
	      q4=q3-ht[i];
	      gradients[p][i]=(-  value(q1, comp)
			       +8*value(q2, comp)
			       -8*value(q3, comp)
			       +  value(q4, comp))/(12*h);
	    }
	break;
      }
      default:
	    Assert(false, ExcInvalidFormula());
    }
}



template <int dim>
void AutoDerivativeFunction<dim>::vector_gradient_list (const typename std::vector<Point<dim> >            &points,
							typename std::vector<std::vector<Tensor<1,dim> > > &gradients) const
{
  Assert (gradients.size() == points.size(),
	  ExcDimensionMismatch(gradients.size(), points.size()));
  for (unsigned p=0; p<points.size(); ++p)
    Assert (gradients[p].size() == n_components,
	    ExcDimensionMismatch(gradients.size(), n_components));
    
  switch (formula)
    {
      case UpwindEuler:
      {
	Tensor<1,dim> q1;
	for (unsigned p=0; p<points.size(); ++p)
	  for (unsigned int i=0; i<dim; ++i)
	    {
	      q1=points[p]-ht[i];
	      for (unsigned int comp=0; comp<n_components; ++comp)
		gradients[p][comp][i]=(value(points[p], comp)-value(q1, comp))/h;
	  }
	break;
      }
      case Euler:
      {
	Tensor<1,dim> q1, q2;
	for (unsigned p=0; p<points.size(); ++p)
	  for (unsigned int i=0; i<dim; ++i)
	    {
	      q1=points[p]+ht[i];
	      q2=points[p]-ht[i];
	      for (unsigned int comp=0; comp<n_components; ++comp)
		gradients[p][comp][i]=(value(q1, comp)-value(q2, comp))/(2*h);
	    }
	break;
      }
      case FourthOrder:
      {
	Tensor<1,dim> q1, q2, q3, q4;
	for (unsigned p=0; p<points.size(); ++p)
	  for (unsigned int i=0; i<dim; ++i)
	    {
	      q2=points[p]+ht[i];
	      q1=q2+ht[i];
	      q3=points[p]-ht[i];
	      q4=q3-ht[i];
	      for (unsigned int comp=0; comp<n_components; ++comp)
		gradients[p][comp][i]=(-  value(q1, comp)
				       +8*value(q2, comp)
				       -8*value(q3, comp)
				       +  value(q4, comp))/(12*h);
	    }
	break;
      }
      default:
	    Assert(false, ExcInvalidFormula());
    }
}


template <int dim>
AutoDerivativeFunction<dim>::DifferenceFormula
AutoDerivativeFunction<dim>::get_formula_of_order(const unsigned int ord)
{
  switch (ord)
    {
      case 0:
      case 1: return UpwindEuler;
      case 2: return Euler;
      case 3:
      case 4: return FourthOrder;
      default:
	    Assert(false, ExcNotImplemented());
    }
  return Euler;
}


template class AutoDerivativeFunction<1>;
template class AutoDerivativeFunction<2>;
template class AutoDerivativeFunction<3>;
