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

#include <math.h>


template <int dim>
TensorProductPolynomials<dim>::TensorProductPolynomials(
  const vector<SmartPointer<Polynomial> > &pols):
		polynomials(pols)
{}

				 // first version of this function
//  template <int dim>
//  void TensorProductPolynomials<dim>::shape_values_and_grads(
//    const Point<dim> &p,
//    vector<double> &values,
//    vector<Tensor<1,dim> > &grads,
//    vector<Tensor<2,dim> > &grad_grads) const
//  {
//    unsigned int n_pols=polynomials->size();
//    unsigned int n_pols2=n_pols*n_pols;
//    unsigned int n_tensor_pols=pow(n_pols, dim);
//    Assert(values.size()==n_tensor_pols || values.size()==0,
//  	 ExcDimensionMismatch2(values.size(), n_tensor_pols, 0));
//    Assert(grads.size()==n_tensor_pols|| grads.size()==0,
//  	 ExcDimensionMismatch2(grads.size(), n_tensor_pols, 0));
//    Assert(grad_grads.size()==n_tensor_pols|| grad_grads.size()==0,
//  	 ExcDimensionMismatch2(grad_grads.size(), n_tensor_pols, 0));

//    unsigned int v_size=0;
//    bool update_values=false, update_grads=false, update_grad_grads=false;
//    if (values.size()==n_tensor_pols)
//      {
//        update_values=true;
//        v_size=1;
//      }
//    if (grads.size()==n_tensor_pols)
//      {
//        update_grads=true;
//        v_size=2;
//      }
//    if (grad_grads.size()==n_tensor_pols)
//      {
//        update_grad_grads=true;
//        v_size=3;
//      }

//    vector<vector<vector<double> > > v(
//      dim, vector<vector<double> > (n_pols, vector<double> (v_size)));

//    for (unsigned int d=0; d<dim; ++d)
//      {
//        vector<vector<double> > &v_d=v[d];
//        Assert(v_d.size()==n_pols, ExcInternalError());
//        for (unsigned int i=0; i<n_pols; ++i)
//  	polynomials->operator[](i)->value(p(d), v_d[i]);
//      }
  
//    if (update_values)
//      {
//  				       // values[i]=p_j(x)
//        for (unsigned i=0; i<n_tensor_pols;)
//  	for (unsigned int j=0; j<n_pols; ++j, ++i)
//  	  values[i]=v[0][j][0];

//  				       // values[i]*=p_j(y)
//        if (dim>1)
//          for (unsigned i=0; i<n_tensor_pols;)
//  	  for (unsigned j=0; j<n_pols; ++j)
//  	    {
//  	      double pjy=v[1][j][0];
//  	      for (unsigned int k=0; k<n_pols; ++k, ++i)
//  		values[i]*=pjy;
//  	    }
      
//  				       // values[i]*=p_j(z)
//        if (dim>2)
//  	for (unsigned int i=0; i<n_tensor_pols;)
//  	  for (unsigned int j=0; j<n_pols; ++j)
//  	    {
//  	      double pjz=v[2][j][0];
//  	      for (unsigned int k=0; k<n_pols2; ++k, ++i)
//  		values[i]*=pjz;
//  	    }
      
//        Assert(dim<=3, ExcNotImplemented());
//      }

      
//    if (update_grads)
//      {
//  				       // grads[i][0]=dp_j(x)
//  				       // grads[i][1]=p_j(x)
//  				       // grads[i][2]=p_j(x)
//        for (unsigned i=0; i<n_tensor_pols;)
//  	for (unsigned int j=0; j<n_pols; ++j, ++i)
//  	  switch (dim)
//  	    {
//  	      case 3: grads[i][2]=v[0][j][0];
//  	      case 2: grads[i][1]=v[0][j][0];
//  	      case 1: grads[i][0]=v[0][j][1]; break;
//  	      default: Assert(false, ExcNotImplemented());
//  	    }


//  				       // grad[i][0]*=p_j(y)
//  				       // grad[i][1]*=dp_j(y)
//  				       // grad[i][2]*=p_j(y)
//        if (dim>1)
//  	for (unsigned i=0; i<n_tensor_pols;)
//  	  for (unsigned j=0; j<n_pols; ++j)
//  	    {
//  	      double pjy=v[1][j][0];
//  	      double dpjy=v[1][j][1];
//  	      for (unsigned int k=0; k<n_pols; ++k, ++i)
//  		{
//  		  grads[i][0]*=pjy;
//  		  grads[i][1]*=dpjy;
//  		  if (dim==3) grads[i][2]*=pjy;
//  		  Assert(dim<=3, ExcNotImplemented());
//  		}
//  	    }
      
//  				       // grad[i][0]*=p_j(z)
//  				       // grad[i][1]*=p_j(z)
//  				       // grad[i][2]*=dp_j(z)
//        if (dim>2)
//  	for (unsigned i=0; i<n_tensor_pols;)
//  	  for (unsigned j=0; j<n_pols; ++j)
//  	    {
//  	      double pjz=v[2][j][0];
//  	      double dpjz=v[2][j][1];
//  	      for (unsigned int k=0; k<n_pols2; ++k, ++i)
//  		{
//  		  grads[i][0]*=pjz;
//  		  grads[i][1]*=pjz;
//  		  grads[i][2]*=dpjz;
//  		  Assert(dim<=3, ExcNotImplemented());
//  		}
//  	    }
//      }

//    if (update_grad_grads)
//      {
//        Assert(false, ExcNotImplemented());
//      }
  
//  }



template <typename number>
number pow(number x, unsigned int y)
{
  number value=1;
  for (unsigned int i=0; i<y; ++i)
    value*=x;
  return value;
}


				 // second and shorter version
template <int dim>
void TensorProductPolynomials<dim>::shape_values_and_grads(
  const Point<dim> &p,
  vector<double> &values,
  vector<Tensor<1,dim> > &grads,
  vector<Tensor<2,dim> > &grad_grads) const
{
  unsigned int n_pols=polynomials.size();
  unsigned int n_tensor_pols=pow(n_pols, dim);
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
  
  vector<unsigned int> n_pols_to(dim+1);
  n_pols_to[0]=1;
  for (unsigned int i=0; i<dim; ++i)
    n_pols_to[i+1]=n_pols_to[i]*n_pols;
  Assert(n_pols_to[dim]==n_tensor_pols, ExcInternalError());
  
  if (update_values)
    {
      for (unsigned int i=0; i<n_tensor_pols; ++i)
	values[i]=1;
      
      for (unsigned int d=0; d<dim; ++d)
	for (unsigned int i=0; i<n_tensor_pols; ++i)
	  values[i]*=v[d][(i/n_pols_to[d])%n_pols][0];
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
