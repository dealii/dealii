//----------------------  tensor_product_polynomials.cc  ------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------  tensor_product_polynomials.cc  ------------


#include <base/exceptions.h>
#include <base/tensor_product_polynomials.h>



//TODO: [RH] Is this function only used inside this class?
// It should be inline, anyway. Also, there is no documentation.
// rather, make it a file-static function, or better in an anonymous
// namespace

template <int dim>
unsigned int TensorProductPolynomials<dim>::power(const unsigned int x,
						  const unsigned int y)
{
  unsigned int value=1;
  for (unsigned int i=0; i<y; ++i)
    value*=x;
  return value;
}



template <int dim>
double
TensorProductPolynomials<dim>::compute_value(const unsigned int i,
					     const Point<dim> &p) const
{
  const unsigned int n_pols=polynomials.size();
  
  double value=1.;
  for (unsigned int d=0; d<dim; ++d)
    value *= polynomials[(i/n_pols_to[d])%n_pols].value(p(d));
  
  return value;
}

  
template <int dim>
Tensor<1,dim>
TensorProductPolynomials<dim>::compute_grad(const unsigned int i,
					    const Point<dim> &p) const
{
  const unsigned int n_pols=polynomials.size();
  
  std::vector<std::vector<double> > v(dim, std::vector<double> (2));

  for (unsigned int d=0; d<dim; ++d)
    polynomials[(i/n_pols_to[d])%n_pols].value(p(d), v[d]);
  
  Tensor<1,dim> grad;
  for (unsigned int d=0; d<dim; ++d)
    grad[d]=1.;
  
  for (unsigned int d=0; d<dim; ++d)
    for (unsigned int x=0; x<dim; ++x)
      grad[d]*=v[x][d==x];
  
  return grad;
}


template <int dim>
Tensor<2,dim>
TensorProductPolynomials<dim>::compute_grad_grad(const unsigned int i,
						 const Point<dim> &p) const
{
  const unsigned int n_pols=polynomials.size();
    
  std::vector<std::vector<double> > v(dim, std::vector<double> (3));
  for (unsigned int d=0; d<dim; ++d)
    polynomials[(i/n_pols_to[d])%n_pols].value(p(d), v[d]);
  
  Tensor<2,dim> grad_grad;

  for (unsigned int d1=0; d1<dim; ++d1)
    for (unsigned int d2=0; d2<dim; ++d2)
      grad_grad[d1][d2]=1.;
  
  for (unsigned int x=0; x<dim; ++x)
    for (unsigned int d1=0; d1<dim; ++d1)
      for (unsigned int d2=0; d2<dim; ++d2)
	{
	  unsigned int derivative=0;
	  if (d1==x || d2==x)
	    {
	      if (d1==d2)
		derivative=2;
	      else
		derivative=1;
	    } 
	  grad_grad[d1][d2]*=v[x][derivative];
	}

  return grad_grad;
}




template <int dim>
void TensorProductPolynomials<dim>::compute(
  const Point<dim>                     &p,
  std::vector<double>                  &values,
  typename std::vector<Tensor<1,dim> > &grads,
  typename std::vector<Tensor<2,dim> > &grad_grads) const
{
  const unsigned int n_pols=polynomials.size();
  
  Assert(values.size()==n_tensor_pols || values.size()==0,
	 ExcDimensionMismatch2(values.size(), n_tensor_pols, 0));
  Assert(grads.size()==n_tensor_pols|| grads.size()==0,
	 ExcDimensionMismatch2(grads.size(), n_tensor_pols, 0));
  Assert(grad_grads.size()==n_tensor_pols|| grad_grads.size()==0,
	 ExcDimensionMismatch2(grad_grads.size(), n_tensor_pols, 0));

  unsigned int v_size=0;
  bool update_values=false, update_grads=false, update_grad_grads=false;
  if (values.size()==n_tensor_pols)
    {
      update_values=true;
      v_size=1;
    }
  if (grads.size()==n_tensor_pols)
    {
      update_grads=true;
      v_size=2;
    }
  if (grad_grads.size()==n_tensor_pols)
    {
      update_grad_grads=true;
      v_size=3;
    }

  std::vector<std::vector<std::vector<double> > > v(
    dim, std::vector<std::vector<double> > (n_pols, std::vector<double> (v_size)));

  for (unsigned int d=0; d<dim; ++d)
    {
      std::vector<std::vector<double> > &v_d=v[d];
      Assert(v_d.size()==n_pols, ExcInternalError());
      for (unsigned int i=0; i<n_pols; ++i)
	polynomials[i].value(p(d), v_d[i]);
    }
  
  if (update_values)
    {
      for (unsigned int i=0; i<n_tensor_pols; ++i)
	values[i]=1;
      
      for (unsigned int x=0; x<dim; ++x)
	for (unsigned int i=0; i<n_tensor_pols; ++i)
	  values[i]*=v[x][(i/n_pols_to[x])%n_pols][0];
    }
  
  if (update_grads)
    {
      for (unsigned int i=0; i<n_tensor_pols; ++i)
	for (unsigned int d=0; d<dim; ++d)
	  grads[i][d]=1.;

      for (unsigned int x=0; x<dim; ++x)
	for (unsigned int i=0; i<n_tensor_pols; ++i)
	  for (unsigned int d=0; d<dim; ++d)
	    grads[i][d]*=v[x][(i/n_pols_to[x])%n_pols][d==x];
    }

  if (update_grad_grads)
    {
      for (unsigned int i=0; i<n_tensor_pols; ++i)
	for (unsigned int d1=0; d1<dim; ++d1)
	  for (unsigned int d2=0; d2<dim; ++d2)
	    grad_grads[i][d1][d2]=1.;

      for (unsigned int x=0; x<dim; ++x)
	for (unsigned int i=0; i<n_tensor_pols; ++i)
	  for (unsigned int d1=0; d1<dim; ++d1)
	    for (unsigned int d2=0; d2<dim; ++d2)
	      {
		unsigned int derivative=0;
		if (d1==x || d2==x)
		  {
		    if (d1==d2)
		      derivative=2;
		    else
		      derivative=1;
		  } 
		grad_grads[i][d1][d2]*=
		  v[x][(i/n_pols_to[x])%n_pols][derivative];
	      }
    }
}


template<int dim> inline
unsigned int
TensorProductPolynomials<dim>::n_tensor_product_polynomials() const
{
  return n_tensor_pols;
}

  
template class TensorProductPolynomials<1>;
template class TensorProductPolynomials<2>;
template class TensorProductPolynomials<3>;
