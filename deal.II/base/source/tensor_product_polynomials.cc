//----------------------  tensor_product_polynomials.cc  ------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------  tensor_product_polynomials.cc  ------------


#include <base/exceptions.h>
#include <base/tensor_product_polynomials.h>



template <int dim>
TensorProductPolynomials<dim>::TensorProductPolynomials(
  const vector<SmartPointer<Polynomial> > &pols):
		polynomials(pols)
{}


template <int dim>
void TensorProductPolynomials<dim>::compute(
  const Point<dim> &p,
  vector<double> &values,
  vector<Tensor<1,dim> > &grads,
  vector<Tensor<2,dim> > &grad_grads) const
{
  unsigned int n_pols=polynomials.size();
  vector<unsigned int> n_pols_to(dim+1);
  n_pols_to[0]=1;
  for (unsigned int i=0; i<dim; ++i)
    n_pols_to[i+1]=n_pols_to[i]*n_pols;
  unsigned int n_tensor_pols=n_pols_to[dim];
  
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

  vector<vector<vector<double> > > v(
    dim, vector<vector<double> > (n_pols, vector<double> (v_size)));

  for (unsigned int d=0; d<dim; ++d)
    {
      vector<vector<double> > &v_d=v[d];
      Assert(v_d.size()==n_pols, ExcInternalError());
      for (unsigned int i=0; i<n_pols; ++i)
	polynomials[i]->value(p(d), v_d[i]);
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



  
template class TensorProductPolynomials<1>;
template class TensorProductPolynomials<2>;
template class TensorProductPolynomials<3>;
