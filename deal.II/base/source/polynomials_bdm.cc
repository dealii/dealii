//-------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------

#include <base/polynomials_bdm.h>
#include <base/quadrature_lib.h>
#include <iostream>
using namespace std;
using namespace Polynomials;


template <int dim>
PolynomialsBDM<dim>::PolynomialsBDM (const unsigned int k)
		:
		polynomial_space (Polynomials::Monomial<double>::generate_complete_basis(k)),
		monomials(1),
		monomial_derivatives(1),
		n_pols(dim * polynomial_space.n()+2),
		p_values(polynomial_space.n()),
		p_grads(polynomial_space.n()),
		p_grad_grads(polynomial_space.n())
{
  Assert (dim == 2, ExcNotImplemented());
  monomials[0] = Monomial<double> (k+1);
  for (unsigned int i=0;i<monomials.size();++i)
    monomial_derivatives[i] = monomials[i].derivative();
}



template <int dim>
void
PolynomialsBDM<dim>::compute (const Point<dim>            &unit_point,
			      std::vector<Tensor<1,dim> > &values,
			      std::vector<Tensor<2,dim> > &grads,
			      std::vector<Tensor<3,dim> > &grad_grads) const
{
  Assert(values.size()==n_pols || values.size()==0,
	 ExcDimensionMismatch2(values.size(), n_pols, 0));
  Assert(grads.size()==n_pols|| grads.size()==0,
	 ExcDimensionMismatch2(grads.size(), n_pols, 0));
  Assert(grad_grads.size()==n_pols|| grad_grads.size()==0,
	 ExcDimensionMismatch2(grad_grads.size(), n_pols, 0));

  const unsigned int n_sub = polynomial_space.n();
  p_values.resize((values.size() == 0) ? 0 : n_sub);
  p_grads.resize((grads.size() == 0) ? 0 : n_sub);
  p_grad_grads.resize((grad_grads.size() == 0) ? 0 : n_sub);

				   // Compute values of complete space
				   // and insert into tensors.  Result
				   // will have first all polynomials
				   // in the x-component, then y and
				   // z.
  polynomial_space.compute (unit_point, p_values, p_grads, p_grad_grads);

  std::fill(values.begin(), values.end(), Tensor<1,dim>());
  for (unsigned int i=0;i<p_values.size();++i)
    {
      for (unsigned int j=0;j<dim;++j)
	{
	  values[i+j*n_sub][j] = p_values[i];
	  std::cerr << i+j*n_sub << ' ' << j << ' ' << p_values[i] << std::endl;
	}
      
    }

				   // Let's hope this is not the transpose
  std::fill(grads.begin(), grads.end(), Tensor<2,dim>());
  for (unsigned int i=0;i<p_grads.size();++i)
    {
      for (unsigned int j=0;j<dim;++j)
	grads[i+j*n_sub][j] = p_grads[i];
    }

				   // Let's hope this is not the transpose
  std::fill(grad_grads.begin(), grad_grads.end(), Tensor<3,dim>());
  for (unsigned int i=0;i<p_grad_grads.size();++i)
    {
      for (unsigned int j=0;j<dim;++j)
	grad_grads[i+j*n_sub][j] = p_grad_grads[i];
    }
  
  const unsigned int start = dim*n_sub;
  if (values.size() != 0)
    {
      values[start][0] = monomials[0].value (unit_point(0));
      values[start][1] = -unit_point(1)
			 * monomial_derivatives[0].value (unit_point(0));
      values[start+1][0] = -unit_point(0)
			 * monomial_derivatives[0].value (unit_point(1));
      values[start+1][1] = monomials[0].value (unit_point(1));
    }
  if (grads.size() != 0)
    {
      Assert(false,ExcNotImplemented());
    }
  if (grad_grads.size() != 0)
    {
      Assert(false,ExcNotImplemented());
    }
}


template <int dim>
void
PolynomialsBDM<dim>::compute_node_matrix (Table<2,double>& A) const
{
  std::vector<Polynomial<double> > legendre(2);
  for (unsigned int i=0;i<legendre.size();++i)
    legendre[i] = Legendre(i);

  QGauss<1> qface(polynomial_space.degree());

  Table<2,double> integrals (n(), n());

  std::vector<Tensor<1,dim> > values(n());
  std::vector<Tensor<2,dim> > grads;
  std::vector<Tensor<3,dim> > grad_grads;
  values.resize(n());

  for (unsigned int face=0;face<2*dim;++face)
    for (unsigned int k=0;k<qface.n_quadrature_points;++k)
      {
	const double w = qface.weight(k);
	const double x = qface.point(k)(0);
	Point<dim> p;
	switch (face)
	  {
	    case 2:
	      p(1) = 1.;
	    case 0:
	      p(0) = x;
	      break;
	    case 1:
	      p(0) = 1.;
	    case 3:
	      p(1) = x;
	      break;	      
	  }
	std::cerr << p << std::endl;
	
	compute (p, values, grads, grad_grads);
	for (unsigned int i=0;i<n();++i)
	  {
	    for (unsigned int j=0;j<legendre.size();++j)
	      A(2*face+j,i) += w * values[i][1-face%2] * legendre[j].value(x);
	  }
      }
				   // Volume integrals are missing
  Assert (polynomial_space.degree() < 2,
	  ExcNotImplemented());
}


template class PolynomialsBDM<2>;

