// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/base/polynomials_abf.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <iostream>
#include <iomanip>


DEAL_II_NAMESPACE_OPEN


template <int dim>
PolynomialsABF<dim>::PolynomialsABF (const unsigned int k)
  :
  my_degree(k),
  polynomial_space (create_polynomials (k)),
  n_pols(compute_n_pols(k))
{
}

template <int dim>
std::vector<std::vector< Polynomials::Polynomial< double > > >
PolynomialsABF<dim>::create_polynomials (const unsigned int k)
{
  std::vector<std::vector< Polynomials::Polynomial< double > > > pols(dim);
  pols[0] = Polynomials::Legendre::generate_complete_basis(k+2);
  if (k == 0)
    for (unsigned int d=1; d<dim; ++d)
      pols[d] = Polynomials::Legendre::generate_complete_basis(0);
  else
    for (unsigned int d=1; d<dim; ++d)
      pols[d] = Polynomials::Legendre::generate_complete_basis(k);

  return pols;
}


template <int dim>
void
PolynomialsABF<dim>::compute (const Point<dim>            &unit_point,
                              std::vector<Tensor<1,dim> > &values,
                              std::vector<Tensor<2,dim> > &grads,
                              std::vector<Tensor<3,dim> > &grad_grads) const
{
  Assert(values.size()==n_pols || values.size()==0,
         ExcDimensionMismatch(values.size(), n_pols));
  Assert(grads.size()==n_pols|| grads.size()==0,
         ExcDimensionMismatch(grads.size(), n_pols));
  Assert(grad_grads.size()==n_pols|| grad_grads.size()==0,
         ExcDimensionMismatch(grad_grads.size(), n_pols));

  // guard access to the scratch
  // arrays in the following block
  // using a mutex to make sure they
  // are not used by multiple threads
  static Threads::Mutex mutex;
  Threads::Mutex::ScopedLock lock(mutex);

  static std::vector<double> p_values;
  static std::vector<Tensor<1,dim> > p_grads;
  static std::vector<Tensor<2,dim> > p_grad_grads;

  const unsigned int n_sub = polynomial_space.n();
  p_values.resize((values.size() == 0) ? 0 : n_sub);
  p_grads.resize((grads.size() == 0) ? 0 : n_sub);
  p_grad_grads.resize((grad_grads.size() == 0) ? 0 : n_sub);

  for (unsigned int d=0; d<dim; ++d)
    {
      // First we copy the point. The
      // polynomial space for
      // component d consists of
      // polynomials of degree k+1 in
      // x_d and degree k in the
      // other variables. in order to
      // simplify this, we use the
      // same AnisotropicPolynomial
      // space and simply rotate the
      // coordinates through all
      // directions.
      Point<dim> p;
      for (unsigned int c=0; c<dim; ++c)
        p(c) = unit_point((c+d)%dim);

      polynomial_space->compute (p, p_values, p_grads, p_grad_grads);

      for (unsigned int i=0; i<p_values.size(); ++i)
        values[i+d*n_sub][d] = p_values[i];

      for (unsigned int i=0; i<p_grads.size(); ++i)
        for (unsigned int d1=0; d1<dim; ++d1)
          grads[i+d*n_sub][d][(d1+d)%dim] = p_grads[i][d1];

      for (unsigned int i=0; i<p_grad_grads.size(); ++i)
        for (unsigned int d1=0; d1<dim; ++d1)
          for (unsigned int d2=0; d2<dim; ++d2)
            grad_grads[i+d*n_sub][d][(d1+d)%dim][(d2+d)%dim]
              = p_grad_grads[i][d1][d2];
    }
}


template <int dim>
unsigned int
PolynomialsABF<dim>::compute_n_pols(unsigned int k)
{
  if (dim == 1) return k+3;
  if (dim == 2) return 2*(k+1)*(k+3);
  //TODO:Check what are the correct numbers ...
  if (dim == 3) return 3*(k+1)*(k+1)*(k+3);
  Assert(false, ExcNotImplemented());
  return 0;
}


template class PolynomialsABF<1>;
template class PolynomialsABF<2>;
template class PolynomialsABF<3>;


DEAL_II_NAMESPACE_CLOSE
