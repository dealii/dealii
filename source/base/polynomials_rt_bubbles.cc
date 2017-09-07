// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


#include "deal.II/base/polynomials_rt_bubbles.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <iostream>
#include <iomanip>

DEAL_II_NAMESPACE_OPEN


template <int dim>
PolynomialsRT_Bubbles<dim>::PolynomialsRT_Bubbles (const unsigned int k)
  :
  my_degree(k),
  polynomial_space (create_polynomials(k-1)),
  monomials(k+2),
  n_pols(compute_n_pols(k))
{
  Assert (dim >= 2, ExcImpossibleInDim(dim));

  for (unsigned int i=0; i<monomials.size(); ++i)
    monomials[i] = Polynomials::Monomial<double> (i);
}


template <int dim>
std::vector<std::vector< Polynomials::Polynomial< double > > >
PolynomialsRT_Bubbles<dim>::create_polynomials (const unsigned int k)
{
  std::vector<std::vector< Polynomials::Polynomial< double > > > pols(dim);
  pols[0] = Polynomials::LagrangeEquidistant::generate_complete_basis(k+1);
  if (k == 0)
    for (unsigned int d=1; d<dim; ++d)
      pols[d] = Polynomials::Legendre::generate_complete_basis(0);
  else
    for (unsigned int d=1; d<dim; ++d)
      pols[d] = Polynomials::LagrangeEquidistant::generate_complete_basis(k);

  return pols;
}

template <int dim>
void
PolynomialsRT_Bubbles<dim>::compute (const Point<dim>         &unit_point,
                                     std::vector<Tensor<1,dim> > &values,
                                     std::vector<Tensor<2,dim> > &grads,
                                     std::vector<Tensor<3,dim> > &grad_grads,
                                     std::vector<Tensor<4,dim> > &third_derivatives,
                                     std::vector<Tensor<5,dim> > &fourth_derivatives) const
{
  Assert(values.size()==n_pols || values.size()==0,
         ExcDimensionMismatch(values.size(), n_pols));
  Assert(grads.size()==n_pols || grads.size()==0,
         ExcDimensionMismatch(grads.size(), n_pols));
  Assert(grad_grads.size()==n_pols || grad_grads.size()==0,
         ExcDimensionMismatch(grad_grads.size(), n_pols));
  Assert(third_derivatives.size()==n_pols || third_derivatives.size()==0,
         ExcDimensionMismatch(third_derivatives.size(), n_pols));
  Assert(fourth_derivatives.size()==n_pols || fourth_derivatives.size()==0,
         ExcDimensionMismatch(fourth_derivatives.size(), n_pols));

  // Third and fourth derivatives are not implemented
  (void)third_derivatives;
  Assert(third_derivatives.size()==0,
         ExcNotImplemented());
  (void)fourth_derivatives;
  Assert(fourth_derivatives.size()==0,
         ExcNotImplemented());

  const unsigned int n_sub = polynomial_space.n();

  // guard access to the scratch arrays in the following block
  // using a mutex to make sure they are not used by multiple threads
  // at once

  static Threads::Mutex mutex;
  Threads::Mutex::ScopedLock lock(mutex);

  static std::vector<double> p_values;
  static std::vector<Tensor<1,dim> > p_grads;
  static std::vector<Tensor<2,dim> > p_grad_grads;
  static std::vector<Tensor<3,dim> > p_third_derivatives;
  static std::vector<Tensor<4,dim> > p_fourth_derivatives;

  p_values.resize((values.size() == 0) ? 0 : n_sub);
  p_grads.resize((grads.size() == 0) ? 0 : n_sub);
  p_grad_grads.resize((grad_grads.size() == 0) ? 0 : n_sub);

  // This is the Raviart-Thomas part of the space
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

      polynomial_space.compute (p, p_values, p_grads, p_grad_grads, p_third_derivatives, p_fourth_derivatives);

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

  // Next we compute the polynomials and derivatives
  // of the curl part of the space
  std::vector<std::vector<double> > monoval_plus(dim, std::vector<double>(4));
  std::vector<std::vector<double> > monoval(dim, std::vector<double>(4));

  std::vector<std::vector<double> > monoval_i(dim, std::vector<double>(4));
  std::vector<std::vector<double> > monoval_j(dim, std::vector<double>(4));
  std::vector<std::vector<double> > monoval_jplus(dim, std::vector<double>(4));

  unsigned int start = dim*n_sub;

  for (unsigned int d=0; d<dim; ++d)
    {
      monomials[my_degree+1].value(unit_point(d), monoval_plus[d]);
      monomials[my_degree].value(unit_point(d), monoval[d]);
    }

  if (dim == 2)
    {
      for (unsigned int i=0; i<=my_degree; ++i, ++start)
        {
          for (unsigned int d=0; d<dim; ++d)
            monomials[i].value(unit_point(d), monoval_i[d]);

          if (values.size() != 0)
            {
              values[start][0] = -monoval_i[0][0] * monoval_plus[1][1];
              values[start][1] = monoval_i[0][1] * monoval_plus[1][0];

              values[start+my_degree+1][0] = monoval_plus[0][0] * monoval_i[1][1];
              values[start+my_degree+1][1] = -monoval_plus[0][1] * monoval_i[1][0];
            }

          if (grads.size() != 0)
            {
              grads[start][0][0] = -monoval_i[0][1] * monoval_plus[1][1];
              grads[start][0][1] = -monoval_i[0][0] * monoval_plus[1][2];
              grads[start][1][0] = monoval_i[0][2] * monoval_plus[1][0];
              grads[start][1][1] = monoval_i[0][1] * monoval_plus[1][1];

              grads[start+my_degree+1][0][0] = monoval_plus[0][1] * monoval_i[1][1];
              grads[start+my_degree+1][0][1] = monoval_plus[0][0] * monoval_i[1][2];
              grads[start+my_degree+1][1][0] = -monoval_plus[0][2] * monoval_i[1][0];
              grads[start+my_degree+1][1][1] = -monoval_plus[0][1] * monoval_i[1][1];
            }

          if (grad_grads.size() != 0)
            {
              grad_grads[start][0][0][0] = -monoval_i[0][2] * monoval_plus[1][1];
              grad_grads[start][0][0][1] = -monoval_i[0][1] * monoval_plus[1][2];
              grad_grads[start][0][1][0] = -monoval_i[0][1] * monoval_plus[1][2];
              grad_grads[start][0][1][1] = -monoval_i[0][0] * monoval_plus[1][3];
              grad_grads[start][1][0][0] = monoval_i[0][3] * monoval_plus[1][0];
              grad_grads[start][1][0][1] = monoval_i[0][2] * monoval_plus[1][1];
              grad_grads[start][1][1][0] = monoval_i[0][2] * monoval_plus[1][1];
              grad_grads[start][1][1][1] = monoval_i[0][1] * monoval_plus[1][2];

              grad_grads[start+my_degree+1][0][0][0] = monoval_plus[0][2] * monoval_i[1][1];
              grad_grads[start+my_degree+1][0][0][1] = monoval_plus[0][1] * monoval_i[1][2];
              grad_grads[start+my_degree+1][0][1][0] = monoval_plus[0][1] * monoval_i[1][2];
              grad_grads[start+my_degree+1][0][1][1] = monoval_plus[0][0] * monoval_i[1][3];
              grad_grads[start+my_degree+1][1][0][0] = -monoval_plus[0][3] * monoval_i[1][0];
              grad_grads[start+my_degree+1][1][0][1] = -monoval_plus[0][2] * monoval_i[1][1];
              grad_grads[start+my_degree+1][1][1][0] = -monoval_plus[0][2] * monoval_i[1][1];
              grad_grads[start+my_degree+1][1][1][1] = -monoval_plus[0][1] * monoval_i[1][2];
            }
        }
      Assert(start == n_pols - my_degree - 1, ExcInternalError());
    }
  else if (dim == 3)
    {
      const unsigned int n_curls = (my_degree+1) * (2*my_degree + 1);
      // Span of $\tilde{B}$
      for (unsigned int i=0; i<=my_degree; ++i)
        {
          for (unsigned int d=0; d<dim; ++d)
            monomials[i].value(unit_point(d), monoval_i[d]);

          for (unsigned int j = 0; j<=my_degree; ++j)
            {
              for (unsigned int d=0; d<dim; ++d)
                {
                  monomials[j].value(unit_point(d), monoval_j[d]);
                  monomials[j+1].value(unit_point(d), monoval_jplus[d]);
                }

              if (values.size() != 0)
                {
                  values[start][0] = monoval_i[0][0] * monoval_j[1][0] * monoval[2][0] * float(j + my_degree + 2.0);
                  values[start][1] = -monoval_i[0][1] * monoval_jplus[1][0] * monoval[2][0];
                  values[start][2] = -monoval_i[0][1] * monoval_j[1][0] * monoval_plus[2][0];

                  values[start+n_curls][0] = -monoval_jplus[0][0] * monoval_i[1][1] * monoval[2][0];
                  values[start+n_curls][1] = monoval_j[0][0] * monoval_i[1][0] * monoval[2][0] * float(j + my_degree + 2.0);
                  values[start+n_curls][2] = -monoval_j[0][0] * monoval_i[1][1] * monoval_plus[2][0];

                  values[start+2*n_curls][0] = -monoval_jplus[0][0] * monoval[1][0] * monoval_i[2][1];
                  values[start+2*n_curls][1] = -monoval_j[0][0] * monoval_plus[1][0] * monoval_i[2][1];
                  values[start+2*n_curls][2] = monoval_j[0][0] * monoval[1][0] * monoval_i[2][0] * float(j + my_degree + 2.0);

                  if (j != my_degree)
                    {
                      values[start+1][0] = monoval_i[0][0] * monoval[1][0] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      values[start+1][1] = -monoval_i[0][1] * monoval_plus[1][0] * monoval_j[2][0];
                      values[start+1][2] = -monoval_i[0][1] * monoval[1][0] * monoval_jplus[2][0];

                      values[start+n_curls+1][0] = -monoval_plus[0][0] * monoval_i[1][1] * monoval_j[2][0];
                      values[start+n_curls+1][1] = monoval[0][0] * monoval_i[1][0] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      values[start+n_curls+1][2] = -monoval[0][0] * monoval_i[1][1] * monoval_jplus[2][0];

                      values[start+2*n_curls+1][0] = -monoval_plus[0][0] * monoval_j[1][0] * monoval_i[2][1];
                      values[start+2*n_curls+1][1] = -monoval[0][0] * monoval_jplus[1][0] * monoval_i[2][1];
                      values[start+2*n_curls+1][2] = monoval[0][0] * monoval_j[1][0] * monoval_i[2][0] * float(j + my_degree + 2.0);
                    }
                }

              if (grads.size() != 0)
                {
                  grads[start][0][0] = monoval_i[0][1] * monoval_j[1][0] * monoval[2][0] * float(j + my_degree + 2.0);
                  grads[start][0][1] = monoval_i[0][0] * monoval_j[1][1] * monoval[2][0] * float(j + my_degree + 2.0);
                  grads[start][0][2] = monoval_i[0][0] * monoval_j[1][0] * monoval[2][1] * float(j + my_degree + 2.0);
                  grads[start][1][0] = -monoval_i[0][2] * monoval_jplus[1][0] * monoval[2][0];
                  grads[start][1][1] = -monoval_i[0][1] * monoval_jplus[1][1] * monoval[2][0];
                  grads[start][1][2] = -monoval_i[0][1] * monoval_jplus[1][0] * monoval[2][1];
                  grads[start][2][0] = -monoval_i[0][2] * monoval_j[1][0] * monoval_plus[2][0];
                  grads[start][2][1] = -monoval_i[0][1] * monoval_j[1][1] * monoval_plus[2][0];
                  grads[start][2][2] = -monoval_i[0][1] * monoval_j[1][0] * monoval_plus[2][1];

                  grads[start+n_curls][0][0] = -monoval_jplus[0][1] * monoval_i[1][1] * monoval[2][0];
                  grads[start+n_curls][0][1] = -monoval_jplus[0][0] * monoval_i[1][2] * monoval[2][0];
                  grads[start+n_curls][0][2] = -monoval_jplus[0][0] * monoval_i[1][1] * monoval[2][1];
                  grads[start+n_curls][1][0] = monoval_j[0][1] * monoval_i[1][0] * monoval[2][0] * float(j + my_degree + 2.0);
                  grads[start+n_curls][1][1] = monoval_j[0][0] * monoval_i[1][1] * monoval[2][0] * float(j + my_degree + 2.0);
                  grads[start+n_curls][1][2] = monoval_j[0][0] * monoval_i[1][0] * monoval[2][1] * float(j + my_degree + 2.0);
                  grads[start+n_curls][2][0] = -monoval_j[0][1] * monoval_i[1][1] * monoval_plus[2][0];
                  grads[start+n_curls][2][1] = -monoval_j[0][0] * monoval_i[1][2] * monoval_plus[2][0];
                  grads[start+n_curls][2][2] = -monoval_j[0][0] * monoval_i[1][1] * monoval_plus[2][1];

                  grads[start+2*n_curls][0][0] = -monoval_jplus[0][1] * monoval[1][0] * monoval_i[2][1];
                  grads[start+2*n_curls][0][1] = -monoval_jplus[0][0] * monoval[1][1] * monoval_i[2][1];
                  grads[start+2*n_curls][0][2] = -monoval_jplus[0][0] * monoval[1][0] * monoval_i[2][2];
                  grads[start+2*n_curls][1][0] = -monoval_j[0][1] * monoval_plus[1][0] * monoval_i[2][1];
                  grads[start+2*n_curls][1][1] = -monoval_j[0][0] * monoval_plus[1][1] * monoval_i[2][1];
                  grads[start+2*n_curls][1][2] = -monoval_j[0][0] * monoval_plus[1][0] * monoval_i[2][2];
                  grads[start+2*n_curls][2][0] = monoval_j[0][1] * monoval[1][0] * monoval_i[2][0] * float(j + my_degree + 2.0);
                  grads[start+2*n_curls][2][1] = monoval_j[0][0] * monoval[1][1] * monoval_i[2][0] * float(j + my_degree + 2.0);
                  grads[start+2*n_curls][2][2] = monoval_j[0][0] * monoval[1][0] * monoval_i[2][1] * float(j + my_degree + 2.0);

                  if (j != my_degree)
                    {
                      grads[start+1][0][0] = monoval_i[0][1] * monoval[1][0] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grads[start+1][0][1] = monoval_i[0][0] * monoval[1][1] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grads[start+1][0][2] = monoval_i[0][0] * monoval[1][0] * monoval_j[2][1] * float(j + my_degree + 2.0);
                      grads[start+1][1][0] = -monoval_i[0][2] * monoval_plus[1][0] * monoval_j[2][0];
                      grads[start+1][1][1] = -monoval_i[0][1] * monoval_plus[1][1] * monoval_j[2][0];
                      grads[start+1][1][2] = -monoval_i[0][1] * monoval_plus[1][0] * monoval_j[2][1];
                      grads[start+1][2][0] = -monoval_i[0][2] * monoval[1][0] * monoval_jplus[2][0];
                      grads[start+1][2][1] = -monoval_i[0][1] * monoval[1][1] * monoval_jplus[2][0];
                      grads[start+1][2][2] = -monoval_i[0][1] * monoval[1][0] * monoval_jplus[2][1];

                      grads[start+n_curls+1][0][0] = -monoval_plus[0][1] * monoval_i[1][1] * monoval_j[2][0];
                      grads[start+n_curls+1][0][1] = -monoval_plus[0][0] * monoval_i[1][2] * monoval_j[2][0];
                      grads[start+n_curls+1][0][2] = -monoval_plus[0][0] * monoval_i[1][1] * monoval_j[2][1];
                      grads[start+n_curls+1][1][0] = monoval[0][1] * monoval_i[1][0] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grads[start+n_curls+1][1][1] = monoval[0][0] * monoval_i[1][1] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grads[start+n_curls+1][1][2] = monoval[0][0] * monoval_i[1][0] * monoval_j[2][1] * float(j + my_degree + 2.0);
                      grads[start+n_curls+1][2][0] = -monoval[0][1] * monoval_i[1][1] * monoval_jplus[2][0];
                      grads[start+n_curls+1][2][1] = -monoval[0][0] * monoval_i[1][2] * monoval_jplus[2][0];
                      grads[start+n_curls+1][2][2] = -monoval[0][0] * monoval_i[1][1] * monoval_jplus[2][1];

                      grads[start+2*n_curls+1][0][0] = -monoval_plus[0][1] * monoval_j[1][0] * monoval_i[2][1];
                      grads[start+2*n_curls+1][0][1] = -monoval_plus[0][0] * monoval_j[1][1] * monoval_i[2][1];
                      grads[start+2*n_curls+1][0][2] = -monoval_plus[0][0] * monoval_j[1][0] * monoval_i[2][2];
                      grads[start+2*n_curls+1][1][0] = -monoval[0][1] * monoval_jplus[1][0] * monoval_i[2][1];
                      grads[start+2*n_curls+1][1][1] = -monoval[0][0] * monoval_jplus[1][1] * monoval_i[2][1];
                      grads[start+2*n_curls+1][1][2] = -monoval[0][0] * monoval_jplus[1][0] * monoval_i[2][2];
                      grads[start+2*n_curls+1][2][0] = monoval[0][1] * monoval_j[1][0] * monoval_i[2][0] * float(j + my_degree + 2.0);
                      grads[start+2*n_curls+1][2][1] = monoval[0][0] * monoval_j[1][1] * monoval_i[2][0] * float(j + my_degree + 2.0);
                      grads[start+2*n_curls+1][2][2] = monoval[0][0] * monoval_j[1][0] * monoval_i[2][1] * float(j + my_degree + 2.0);
                    }
                }

              if (grad_grads.size() != 0)
                {
                  grad_grads[start][0][0][0] = monoval_i[0][2] * monoval_j[1][0] * monoval[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start][0][0][1] = monoval_i[0][1] * monoval_j[1][1] * monoval[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start][0][0][2] = monoval_i[0][1] * monoval_j[1][0] * monoval[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start][0][1][0] = monoval_i[0][1] * monoval_j[1][1] * monoval[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start][0][1][1] = monoval_i[0][0] * monoval_j[1][2] * monoval[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start][0][1][2] = monoval_i[0][0] * monoval_j[1][1] * monoval[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start][0][2][0] = monoval_i[0][1] * monoval_j[1][0] * monoval[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start][0][2][1] = monoval_i[0][0] * monoval_j[1][1] * monoval[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start][0][2][2] = monoval_i[0][0] * monoval_j[1][0] * monoval[2][2] * float(j + my_degree + 2.0);
                  grad_grads[start][1][0][0] = -monoval_i[0][3] * monoval_jplus[1][0] * monoval[2][0];
                  grad_grads[start][1][0][1] = -monoval_i[0][2] * monoval_jplus[1][1] * monoval[2][0];
                  grad_grads[start][1][0][2] = -monoval_i[0][2] * monoval_jplus[1][0] * monoval[2][1];
                  grad_grads[start][1][1][0] = -monoval_i[0][2] * monoval_jplus[1][1] * monoval[2][0];
                  grad_grads[start][1][1][1] = -monoval_i[0][1] * monoval_jplus[1][2] * monoval[2][0];
                  grad_grads[start][1][1][2] = -monoval_i[0][1] * monoval_jplus[1][1] * monoval[2][1];
                  grad_grads[start][1][2][0] = -monoval_i[0][2] * monoval_jplus[1][0] * monoval[2][1];
                  grad_grads[start][1][2][1] = -monoval_i[0][1] * monoval_jplus[1][1] * monoval[2][1];
                  grad_grads[start][1][2][2] = -monoval_i[0][1] * monoval_jplus[1][0] * monoval[2][2];
                  grad_grads[start][2][0][0] = -monoval_i[0][3] * monoval_j[1][0] * monoval_plus[2][0];
                  grad_grads[start][2][0][1] = -monoval_i[0][2] * monoval_j[1][1] * monoval_plus[2][0];
                  grad_grads[start][2][0][2] = -monoval_i[0][2] * monoval_j[1][0] * monoval_plus[2][1];
                  grad_grads[start][2][1][0] = -monoval_i[0][2] * monoval_j[1][1] * monoval_plus[2][0];
                  grad_grads[start][2][1][1] = -monoval_i[0][1] * monoval_j[1][2] * monoval_plus[2][0];
                  grad_grads[start][2][1][2] = -monoval_i[0][1] * monoval_j[1][1] * monoval_plus[2][1];
                  grad_grads[start][2][2][0] = -monoval_i[0][2] * monoval_j[1][0] * monoval_plus[2][1];
                  grad_grads[start][2][2][1] = -monoval_i[0][1] * monoval_j[1][1] * monoval_plus[2][1];
                  grad_grads[start][2][2][2] = -monoval_i[0][1] * monoval_j[1][0] * monoval_plus[2][2];

                  grad_grads[start+n_curls][0][0][0] = -monoval_jplus[0][2] * monoval_i[1][1] * monoval[2][0];
                  grad_grads[start+n_curls][0][0][1] = -monoval_jplus[0][1] * monoval_i[1][2] * monoval[2][0];
                  grad_grads[start+n_curls][0][0][2] = -monoval_jplus[0][1] * monoval_i[1][1] * monoval[2][1];
                  grad_grads[start+n_curls][0][1][0] = -monoval_jplus[0][1] * monoval_i[1][2] * monoval[2][0];
                  grad_grads[start+n_curls][0][1][1] = -monoval_jplus[0][0] * monoval_i[1][3] * monoval[2][0];
                  grad_grads[start+n_curls][0][1][2] = -monoval_jplus[0][0] * monoval_i[1][2] * monoval[2][1];
                  grad_grads[start+n_curls][0][2][0] = -monoval_jplus[0][1] * monoval_i[1][1] * monoval[2][1];
                  grad_grads[start+n_curls][0][2][1] = -monoval_jplus[0][0] * monoval_i[1][2] * monoval[2][1];
                  grad_grads[start+n_curls][0][2][2] = -monoval_jplus[0][0] * monoval_i[1][1] * monoval[2][2];
                  grad_grads[start+n_curls][1][0][0] = monoval_j[0][2] * monoval_i[1][0] * monoval[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start+n_curls][1][0][1] = monoval_j[0][1] * monoval_i[1][1] * monoval[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start+n_curls][1][0][2] = monoval_j[0][1] * monoval_i[1][0] * monoval[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start+n_curls][1][1][0] = monoval_j[0][1] * monoval_i[1][1] * monoval[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start+n_curls][1][1][1] = monoval_j[0][0] * monoval_i[1][2] * monoval[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start+n_curls][1][1][2] = monoval_j[0][0] * monoval_i[1][1] * monoval[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start+n_curls][1][2][0] = monoval_j[0][1] * monoval_i[1][0] * monoval[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start+n_curls][1][2][1] = monoval_j[0][0] * monoval_i[1][1] * monoval[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start+n_curls][1][2][2] = monoval_j[0][0] * monoval_i[1][0] * monoval[2][2] * float(j + my_degree + 2.0);
                  grad_grads[start+n_curls][2][0][0] = -monoval_j[0][2] * monoval_i[1][1] * monoval_plus[2][0];
                  grad_grads[start+n_curls][2][0][1] = -monoval_j[0][1] * monoval_i[1][2] * monoval_plus[2][0];
                  grad_grads[start+n_curls][2][0][2] = -monoval_j[0][1] * monoval_i[1][1] * monoval_plus[2][1];
                  grad_grads[start+n_curls][2][1][0] = -monoval_j[0][1] * monoval_i[1][2] * monoval_plus[2][0];
                  grad_grads[start+n_curls][2][1][1] = -monoval_j[0][0] * monoval_i[1][3] * monoval_plus[2][0];
                  grad_grads[start+n_curls][2][1][2] = -monoval_j[0][0] * monoval_i[1][2] * monoval_plus[2][1];
                  grad_grads[start+n_curls][2][2][0] = -monoval_j[0][1] * monoval_i[1][1] * monoval_plus[2][1];
                  grad_grads[start+n_curls][2][2][1] = -monoval_j[0][0] * monoval_i[1][2] * monoval_plus[2][1];
                  grad_grads[start+n_curls][2][2][2] = -monoval_j[0][0] * monoval_i[1][1] * monoval_plus[2][2];

                  grad_grads[start+2*n_curls][0][0][0] = -monoval_jplus[0][2] * monoval[1][0] * monoval_i[2][1];
                  grad_grads[start+2*n_curls][0][0][1] = -monoval_jplus[0][1] * monoval[1][1] * monoval_i[2][1];
                  grad_grads[start+2*n_curls][0][0][2] = -monoval_jplus[0][1] * monoval[1][0] * monoval_i[2][2];
                  grad_grads[start+2*n_curls][0][1][0] = -monoval_jplus[0][1] * monoval[1][1] * monoval_i[2][1];
                  grad_grads[start+2*n_curls][0][1][1] = -monoval_jplus[0][0] * monoval[1][2] * monoval_i[2][1];
                  grad_grads[start+2*n_curls][0][1][2] = -monoval_jplus[0][0] * monoval[1][1] * monoval_i[2][2];
                  grad_grads[start+2*n_curls][0][2][0] = -monoval_jplus[0][1] * monoval[1][0] * monoval_i[2][2];
                  grad_grads[start+2*n_curls][0][2][1] = -monoval_jplus[0][0] * monoval[1][1] * monoval_i[2][2];
                  grad_grads[start+2*n_curls][0][2][2] = -monoval_jplus[0][0] * monoval[1][0] * monoval_i[2][3];
                  grad_grads[start+2*n_curls][1][0][0] = -monoval_j[0][2] * monoval_plus[1][0] * monoval_i[2][1];
                  grad_grads[start+2*n_curls][1][0][1] = -monoval_j[0][1] * monoval_plus[1][1] * monoval_i[2][1];
                  grad_grads[start+2*n_curls][1][0][2] = -monoval_j[0][1] * monoval_plus[1][0] * monoval_i[2][2];
                  grad_grads[start+2*n_curls][1][1][0] = -monoval_j[0][1] * monoval_plus[1][1] * monoval_i[2][1];
                  grad_grads[start+2*n_curls][1][1][1] = -monoval_j[0][0] * monoval_plus[1][2] * monoval_i[2][1];
                  grad_grads[start+2*n_curls][1][1][2] = -monoval_j[0][0] * monoval_plus[1][1] * monoval_i[2][2];
                  grad_grads[start+2*n_curls][1][2][0] = -monoval_j[0][1] * monoval_plus[1][0] * monoval_i[2][2];
                  grad_grads[start+2*n_curls][1][2][1] = -monoval_j[0][0] * monoval_plus[1][1] * monoval_i[2][2];
                  grad_grads[start+2*n_curls][1][2][2] = -monoval_j[0][0] * monoval_plus[1][0] * monoval_i[2][3];
                  grad_grads[start+2*n_curls][2][0][0] = monoval_j[0][2] * monoval[1][0] * monoval_i[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start+2*n_curls][2][0][1] = monoval_j[0][1] * monoval[1][1] * monoval_i[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start+2*n_curls][2][0][2] = monoval_j[0][1] * monoval[1][0] * monoval_i[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start+2*n_curls][2][1][0] = monoval_j[0][1] * monoval[1][1] * monoval_i[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start+2*n_curls][2][1][1] = monoval_j[0][0] * monoval[1][2] * monoval_i[2][0] * float(j + my_degree + 2.0);
                  grad_grads[start+2*n_curls][2][1][2] = monoval_j[0][0] * monoval[1][1] * monoval_i[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start+2*n_curls][2][2][0] = monoval_j[0][1] * monoval[1][0] * monoval_i[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start+2*n_curls][2][2][1] = monoval_j[0][0] * monoval[1][1] * monoval_i[2][1] * float(j + my_degree + 2.0);
                  grad_grads[start+2*n_curls][2][2][2] = monoval_j[0][0] * monoval[1][0] * monoval_i[2][2] * float(j + my_degree + 2.0);

                  if (j != my_degree)
                    {
                      grad_grads[start+1][0][0][0] = monoval_i[0][2] * monoval[1][0] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+1][0][0][1] = monoval_i[0][1] * monoval[1][1] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+1][0][0][2] = monoval_i[0][1] * monoval[1][0] * monoval_j[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+1][0][1][0] = monoval_i[0][1] * monoval[1][1] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+1][0][1][1] = monoval_i[0][0] * monoval[1][2] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+1][0][1][2] = monoval_i[0][0] * monoval[1][1] * monoval_j[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+1][0][2][0] = monoval_i[0][1] * monoval[1][0] * monoval_j[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+1][0][2][1] = monoval_i[0][0] * monoval[1][1] * monoval_j[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+1][0][2][2] = monoval_i[0][0] * monoval[1][0] * monoval_j[2][2] * float(j + my_degree + 2.0);
                      grad_grads[start+1][1][0][0] = -monoval_i[0][3] * monoval_plus[1][0] * monoval_j[2][0];
                      grad_grads[start+1][1][0][1] = -monoval_i[0][2] * monoval_plus[1][1] * monoval_j[2][0];
                      grad_grads[start+1][1][0][2] = -monoval_i[0][2] * monoval_plus[1][0] * monoval_j[2][1];
                      grad_grads[start+1][1][1][0] = -monoval_i[0][2] * monoval_plus[1][1] * monoval_j[2][0];
                      grad_grads[start+1][1][1][1] = -monoval_i[0][1] * monoval_plus[1][2] * monoval_j[2][0];
                      grad_grads[start+1][1][1][2] = -monoval_i[0][1] * monoval_plus[1][1] * monoval_j[2][1];
                      grad_grads[start+1][1][2][0] = -monoval_i[0][2] * monoval_plus[1][0] * monoval_j[2][1];
                      grad_grads[start+1][1][2][1] = -monoval_i[0][1] * monoval_plus[1][1] * monoval_j[2][1];
                      grad_grads[start+1][1][2][2] = -monoval_i[0][1] * monoval_plus[1][0] * monoval_j[2][2];
                      grad_grads[start+1][2][0][0] = -monoval_i[0][3] * monoval[1][0] * monoval_jplus[2][0];
                      grad_grads[start+1][2][0][1] = -monoval_i[0][2] * monoval[1][1] * monoval_jplus[2][0];
                      grad_grads[start+1][2][0][2] = -monoval_i[0][2] * monoval[1][0] * monoval_jplus[2][1];
                      grad_grads[start+1][2][1][0] = -monoval_i[0][2] * monoval[1][1] * monoval_jplus[2][0];
                      grad_grads[start+1][2][1][1] = -monoval_i[0][1] * monoval[1][2] * monoval_jplus[2][0];
                      grad_grads[start+1][2][1][2] = -monoval_i[0][1] * monoval[1][1] * monoval_jplus[2][1];
                      grad_grads[start+1][2][2][0] = -monoval_i[0][2] * monoval[1][0] * monoval_jplus[2][1];
                      grad_grads[start+1][2][2][1] = -monoval_i[0][1] * monoval[1][1] * monoval_jplus[2][1];
                      grad_grads[start+1][2][2][2] = -monoval_i[0][1] * monoval[1][0] * monoval_jplus[2][2];

                      grad_grads[start+n_curls+1][0][0][0] = -monoval_plus[0][2] * monoval_i[1][1] * monoval_j[2][0];
                      grad_grads[start+n_curls+1][0][0][1] = -monoval_plus[0][1] * monoval_i[1][2] * monoval_j[2][0];
                      grad_grads[start+n_curls+1][0][0][2] = -monoval_plus[0][1] * monoval_i[1][1] * monoval_j[2][1];
                      grad_grads[start+n_curls+1][0][1][0] = -monoval_plus[0][1] * monoval_i[1][2] * monoval_j[2][0];
                      grad_grads[start+n_curls+1][0][1][1] = -monoval_plus[0][0] * monoval_i[1][3] * monoval_j[2][0];
                      grad_grads[start+n_curls+1][0][1][2] = -monoval_plus[0][0] * monoval_i[1][2] * monoval_j[2][1];
                      grad_grads[start+n_curls+1][0][2][0] = -monoval_plus[0][1] * monoval_i[1][1] * monoval_j[2][1];
                      grad_grads[start+n_curls+1][0][2][1] = -monoval_plus[0][0] * monoval_i[1][2] * monoval_j[2][1];
                      grad_grads[start+n_curls+1][0][2][2] = -monoval_plus[0][0] * monoval_i[1][1] * monoval_j[2][2];
                      grad_grads[start+n_curls+1][1][0][0] = monoval[0][2] * monoval_i[1][0] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+n_curls+1][1][0][1] = monoval[0][1] * monoval_i[1][1] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+n_curls+1][1][0][2] = monoval[0][1] * monoval_i[1][0] * monoval_j[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+n_curls+1][1][1][0] = monoval[0][1] * monoval_i[1][1] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+n_curls+1][1][1][1] = monoval[0][0] * monoval_i[1][2] * monoval_j[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+n_curls+1][1][1][2] = monoval[0][0] * monoval_i[1][1] * monoval_j[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+n_curls+1][1][2][0] = monoval[0][1] * monoval_i[1][0] * monoval_j[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+n_curls+1][1][2][1] = monoval[0][0] * monoval_i[1][1] * monoval_j[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+n_curls+1][1][2][2] = monoval[0][0] * monoval_i[1][0] * monoval_j[2][2] * float(j + my_degree + 2.0);
                      grad_grads[start+n_curls+1][2][0][0] = -monoval[0][2] * monoval_i[1][1] * monoval_jplus[2][0];
                      grad_grads[start+n_curls+1][2][0][1] = -monoval[0][1] * monoval_i[1][2] * monoval_jplus[2][0];
                      grad_grads[start+n_curls+1][2][0][2] = -monoval[0][1] * monoval_i[1][1] * monoval_jplus[2][1];
                      grad_grads[start+n_curls+1][2][1][0] = -monoval[0][1] * monoval_i[1][2] * monoval_jplus[2][0];
                      grad_grads[start+n_curls+1][2][1][1] = -monoval[0][0] * monoval_i[1][3] * monoval_jplus[2][0];
                      grad_grads[start+n_curls+1][2][1][2] = -monoval[0][0] * monoval_i[1][2] * monoval_jplus[2][1];
                      grad_grads[start+n_curls+1][2][2][0] = -monoval[0][1] * monoval_i[1][1] * monoval_jplus[2][1];
                      grad_grads[start+n_curls+1][2][2][1] = -monoval[0][0] * monoval_i[1][2] * monoval_jplus[2][1];
                      grad_grads[start+n_curls+1][2][2][2] = -monoval[0][0] * monoval_i[1][1] * monoval_jplus[2][2];

                      grad_grads[start+2*n_curls+1][0][0][0] = -monoval_plus[0][2] * monoval_j[1][0] * monoval_i[2][1];
                      grad_grads[start+2*n_curls+1][0][0][1] = -monoval_plus[0][1] * monoval_j[1][1] * monoval_i[2][1];
                      grad_grads[start+2*n_curls+1][0][0][2] = -monoval_plus[0][1] * monoval_j[1][0] * monoval_i[2][2];
                      grad_grads[start+2*n_curls+1][0][1][0] = -monoval_plus[0][1] * monoval_j[1][1] * monoval_i[2][1];
                      grad_grads[start+2*n_curls+1][0][1][1] = -monoval_plus[0][0] * monoval_j[1][2] * monoval_i[2][1];
                      grad_grads[start+2*n_curls+1][0][1][2] = -monoval_plus[0][0] * monoval_j[1][1] * monoval_i[2][2];
                      grad_grads[start+2*n_curls+1][0][2][0] = -monoval_plus[0][1] * monoval_j[1][0] * monoval_i[2][2];
                      grad_grads[start+2*n_curls+1][0][2][1] = -monoval_plus[0][0] * monoval_j[1][1] * monoval_i[2][2];
                      grad_grads[start+2*n_curls+1][0][2][2] = -monoval_plus[0][0] * monoval_j[1][0] * monoval_i[2][3];
                      grad_grads[start+2*n_curls+1][1][0][0] = -monoval[0][2] * monoval_jplus[1][0] * monoval_i[2][1];
                      grad_grads[start+2*n_curls+1][1][0][1] = -monoval[0][1] * monoval_jplus[1][1] * monoval_i[2][1];
                      grad_grads[start+2*n_curls+1][1][0][2] = -monoval[0][1] * monoval_jplus[1][0] * monoval_i[2][2];
                      grad_grads[start+2*n_curls+1][1][1][0] = -monoval[0][1] * monoval_jplus[1][1] * monoval_i[2][1];
                      grad_grads[start+2*n_curls+1][1][1][1] = -monoval[0][0] * monoval_jplus[1][2] * monoval_i[2][1];
                      grad_grads[start+2*n_curls+1][1][1][2] = -monoval[0][0] * monoval_jplus[1][1] * monoval_i[2][2];
                      grad_grads[start+2*n_curls+1][1][2][0] = -monoval[0][1] * monoval_jplus[1][0] * monoval_i[2][2];
                      grad_grads[start+2*n_curls+1][1][2][1] = -monoval[0][0] * monoval_jplus[1][1] * monoval_i[2][2];
                      grad_grads[start+2*n_curls+1][1][2][2] = -monoval[0][0] * monoval_jplus[1][0] * monoval_i[2][3];
                      grad_grads[start+2*n_curls+1][2][0][0] = monoval[0][2] * monoval_j[1][0] * monoval_i[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+2*n_curls+1][2][0][1] = monoval[0][1] * monoval_j[1][1] * monoval_i[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+2*n_curls+1][2][0][2] = monoval[0][1] * monoval_j[1][0] * monoval_i[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+2*n_curls+1][2][1][0] = monoval[0][1] * monoval_j[1][1] * monoval_i[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+2*n_curls+1][2][1][1] = monoval[0][0] * monoval_j[1][2] * monoval_i[2][0] * float(j + my_degree + 2.0);
                      grad_grads[start+2*n_curls+1][2][1][2] = monoval[0][0] * monoval_j[1][1] * monoval_i[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+2*n_curls+1][2][2][0] = monoval[0][1] * monoval_j[1][0] * monoval_i[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+2*n_curls+1][2][2][1] = monoval[0][0] * monoval_j[1][1] * monoval_i[2][1] * float(j + my_degree + 2.0);
                      grad_grads[start+2*n_curls+1][2][2][2] = monoval[0][0] * monoval_j[1][0] * monoval_i[2][2] * float(j + my_degree + 2.0);
                    }
                }

              if (j == my_degree)
                start += 1;
              else
                start += 2;

            }
        }
      Assert(start == n_pols - 2*n_curls, ExcInternalError());
    }
}


template <int dim>
unsigned int
PolynomialsRT_Bubbles<dim>::compute_n_pols(unsigned int k)
{
  if (dim == 1 || dim == 2 || dim == 3)
    return dim * pow(k+1, dim);

  Assert(false, ExcNotImplemented());
  return 0;
}

template class PolynomialsRT_Bubbles<1>;
template class PolynomialsRT_Bubbles<2>;
template class PolynomialsRT_Bubbles<3>;


DEAL_II_NAMESPACE_CLOSE
