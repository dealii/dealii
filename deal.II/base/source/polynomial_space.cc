//----------------------  polynomial_space.cc  ------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------  polynomial_space.cc  ------------


#include <base/exceptions.h>
#include <base/polynomial_space.h>



template <int dim>
unsigned int PolynomialSpace<dim>::power(const unsigned int x,
				     const unsigned int y)
{
  unsigned int value=1;
  for (unsigned int i=0; i<y; ++i)
    value*=x;
  return value;
}



template <int dim>
double
PolynomialSpace<dim>::compute_value(const unsigned int /*i*/,
				    const Point<dim> & /*p*/) const
{
  Assert(false, ExcNotImplemented());
  return 0.;
}

  
template <int dim>
Tensor<1,dim>
PolynomialSpace<dim>::compute_grad(const unsigned int /*i*/,
				   const Point<dim> &/*p*/) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1,dim>();
}


template <int dim>
Tensor<2,dim>
PolynomialSpace<dim>::compute_grad_grad(const unsigned int /*i*/,
					const Point<dim> &/*p*/) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<2,dim>();
}




template <int dim>
void PolynomialSpace<dim>::compute(
  const Point<dim>                     &p,
  std::vector<double>                  &values,
  typename std::vector<Tensor<1,dim> > &grads,
  typename std::vector<Tensor<2,dim> > &grad_grads) const
{
  const unsigned int n_1d=polynomials.size();
  
  Assert(values.size()==n_pols || values.size()==0,
	 ExcDimensionMismatch2(values.size(), n_pols, 0));
  Assert(grads.size()==n_pols|| grads.size()==0,
	 ExcDimensionMismatch2(grads.size(), n_pols, 0));
  Assert(grad_grads.size()==n_pols|| grad_grads.size()==0,
	 ExcDimensionMismatch2(grad_grads.size(), n_pols, 0));

  unsigned int v_size=0;
  bool update_values=false, update_grads=false, update_grad_grads=false;
  if (values.size()==n_pols)
    {
      update_values=true;
      v_size=1;
    }
  if (grads.size()==n_pols)
    {
      update_grads=true;
      v_size=2;
    }
  if (grad_grads.size()==n_pols)
    {
      update_grad_grads=true;
      v_size=3;
    }

				   // Store data in a single
				   // vector. Access is by
				   // v[d][n][o]
				   //  d: coordinate direction
				   //  n: number of 1d polynomial
				   //  o: order of derivative
  std::vector<std::vector<std::vector<double> > >
    v(dim,
      std::vector<std::vector<double> > (n_1d,
					 std::vector<double> (v_size, 0.)));

  for (unsigned int d=0; d<dim; ++d)
    {
      std::vector<std::vector<double> >& v_d=v[d];
      Assert(v_d.size()==n_1d, ExcInternalError());
      for (unsigned int i=0; i<n_1d; ++i)
	polynomials[i].value(p(d), v_d[i]);
    }

  if (update_values)
    {
      unsigned int k = 0;
      
      for (unsigned int iz=0;iz<((dim>2) ? n_1d : 1);++iz)
	for (unsigned int iy=0;iy<((dim>1) ? n_1d-iz : 1);++iy)
	  for (unsigned int ix=0; ix<n_1d-iy-iz; ++ix)
	    values[k++] = v[0][ix][0]
			  * ((dim>1) ? v[1][iy][0] : 1.)
			  * ((dim>2) ? v[2][iz][0] : 1.);
    }
  
  if (update_grads)
    {
      unsigned int k = 0;
      
      for (unsigned int iz=0;iz<((dim>2) ? n_1d : 1);++iz)
	for (unsigned int iy=0;iy<((dim>1) ? n_1d-iz : 1);++iy)
	  for (unsigned int ix=0; ix<n_1d-iy-iz; ++ix)
	    {
	      for (unsigned int d=0;d<dim;++d)
		grads[k][d] = v[0][ix][(d==0) ? 1 : 0]
    * ((dim>1) ? v[1][iy][(d==1) ? 1 : 0] : 1.)
    * ((dim>2) ? v[2][iz][(d==2) ? 1 : 0] : 1.);
	      ++k;
	    }
    }

  if (update_grad_grads)
    {
      unsigned int k = 0;
      
      for (unsigned int iz=0;iz<((dim>2) ? n_1d : 1);++iz)
	for (unsigned int iy=0;iy<((dim>1) ? n_1d-iz : 1);++iy)
	  for (unsigned int ix=0; ix<n_1d-iy-iz; ++ix)
	    {
	      for (unsigned int d1=0; d1<dim; ++d1)
		for (unsigned int d2=0; d2<dim; ++d2)
		  {
						     // Derivative
						     // order for each
						     // direction
		    const unsigned int
		      j0 = ((d1==0) ? 1 : 0) + ((d2==0) ? 1 : 0);
		    const unsigned int
		      j1 = ((d1==1) ? 1 : 0) + ((d2==1) ? 1 : 0);
		    const unsigned int
		      j2 = ((d1==2) ? 1 : 0) + ((d2==2) ? 1 : 0);
		    
		    grad_grads[k][d1][d2] = v[0][ix][j0]
					    * ((dim>1) ? v[1][iy][j1] : 1.)
					    * ((dim>2) ? v[2][iz][j2] : 1.);
		  }
	      ++k;
	    }
    }
}


template<int dim>
unsigned int
PolynomialSpace<dim>::n() const
{
  return n_pols;
}

  
template class PolynomialSpace<1>;
template class PolynomialSpace<2>;
template class PolynomialSpace<3>;
