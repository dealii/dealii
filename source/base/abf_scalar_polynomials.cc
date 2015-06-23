// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
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

#include <deal.II/base/abf_scalar_polynomials.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>

DEAL_II_NAMESPACE_OPEN



/* ------------------- ABFScalarPolynomials -------------- */


namespace internal
{
  namespace
  {
    template <int dim>
    inline
    void compute_tensor_index(const unsigned int,
                              const unsigned int,
                              const unsigned int,
                              unsigned int      ( &)[dim])
    {
      Assert(false, ExcNotImplemented());
    }

    inline
    void compute_tensor_index(const unsigned int n,
                              const unsigned int ,
                              const unsigned int ,
                              unsigned int       (&indices)[1])
    {
      indices[0] = n;
    }

    inline
    void compute_tensor_index(const unsigned int n,
                              const unsigned int n_pols_0,
                              const unsigned int ,
                              unsigned int       (&indices)[2])
    {
      if(n < n_pols_0*n_pols_0)
      {
        indices[0] = n % n_pols_0;
        indices[1] = n / n_pols_0;
      }else{
        if(n < n_pols_0*n_pols_0+n_pols_0){
          indices[0] = n_pols_0;
          indices[1] = n % n_pols_0;
        }else{
          indices[0] = n % n_pols_0;
          indices[1] = n_pols_0;
        }
      }
    }

    inline
    void compute_tensor_index(const unsigned int n,
                              const unsigned int n_pols_0,
                              const unsigned int n_pols_1,
                              unsigned int       (&indices)[3])
    {
      if(n < n_pols_0 * n_pols_0 *n_pols_0)
      {
        indices[0] = n % n_pols_0;
        indices[1] = (n/n_pols_0) % n_pols_1;
        indices[2] = n / (n_pols_0*n_pols_1);
      }else{
        if(n < n_pols_0 * n_pols_0 *n_pols_0 + n_pols_0 * n_pols_0)
        {
          indices[0] = n_pols_0;
          indices[1] = n % n_pols_0;
          indices[2] = (n/n_pols_0) % n_pols_1;
        }else{
          if(n < n_pols_0 * n_pols_0 *n_pols_0 + 2*(n_pols_0 * n_pols_0))
          {
            indices[0] = n % n_pols_0;
            indices[1] = n_pols_0;
            indices[2] = (n/n_pols_0) % n_pols_1;
          }else{
            indices[0] = n % n_pols_0;
            indices[1] = (n/n_pols_0) % n_pols_1;
            indices[2] = n_pols_0;
          }
        }
      }
      
    }
  }
}



template <int dim, typename POLY>
inline
void
ABFScalarPolynomials<dim,POLY>::
compute_index (const unsigned int i,
               unsigned int       (&indices)[(dim > 0 ? dim : 1)]) const
{
  Assert (i<n_tensor_pols,ExcInternalError());
  internal::compute_tensor_index(index_map[i], polynomials.size()-1,
                                 polynomials.size()-1, indices);
}



template <int dim, typename POLY>
void
ABFScalarPolynomials<dim,POLY>::output_indices(std::ostream &out) const
{
  unsigned int ix[dim];
  for (unsigned int i=0; i<n_tensor_pols; ++i)
    {
      compute_index(i,ix);
      out << i << "\t";
      for (unsigned int d=0; d<dim; ++d)
        out << ix[d] << " ";
      out << std::endl;
    }
}



template <int dim, typename POLY>
void
ABFScalarPolynomials<dim,POLY>::set_numbering(
  const std::vector<unsigned int> &renumber)
{
  Assert(renumber.size()==index_map.size(),
         ExcDimensionMismatch(renumber.size(), index_map.size()));

  index_map=renumber;
  for (unsigned int i=0; i<index_map.size(); ++i)
    index_map_inverse[index_map[i]]=i;
}



template <>
double
ABFScalarPolynomials<0,Polynomials::Polynomial<double> >
::compute_value(const unsigned int,
                const Point<0> &) const
{
  Assert (false, ExcNotImplemented());
  return 0;
}



template <int dim, typename POLY>
double
ABFScalarPolynomials<dim,POLY>::compute_value (const unsigned int i,
                                                   const Point<dim> &p) const
{
  Assert(dim>0, ExcNotImplemented());

  unsigned int indices[dim];
  compute_index (i, indices);

  double value=1.;
  for (unsigned int d=0; d<dim; ++d)
    value *= polynomials[indices[d]].value(p(d));

  return value;
}



template <int dim, typename POLY>
Tensor<1,dim>
ABFScalarPolynomials<dim,POLY>::compute_grad (const unsigned int i,
                                                  const Point<dim> &p) const
{
  unsigned int indices[dim];
  compute_index (i, indices);

  // compute values and
  // uni-directional derivatives at
  // the given point in each
  // co-ordinate direction
  double v [dim][2];
  {
    std::vector<double> tmp (2);
    for (unsigned int d=0; d<dim; ++d)
      {
        polynomials[indices[d]].value (p(d), tmp);
        v[d][0] = tmp[0];
        v[d][1] = tmp[1];
      }
  }

  Tensor<1,dim> grad;
  for (unsigned int d=0; d<dim; ++d)
    {
      grad[d] = 1.;
      for (unsigned int x=0; x<dim; ++x)
        grad[d] *= v[x][d==x];
    }

  return grad;
}



template <int dim, typename POLY>
Tensor<2,dim>
ABFScalarPolynomials<dim,POLY>::compute_grad_grad (const unsigned int i,
                                                       const Point<dim> &p) const
{
  unsigned int indices[dim];
  compute_index (i, indices);

  double v [dim][3];
  {
    std::vector<double> tmp (3);
    for (unsigned int d=0; d<dim; ++d)
      {
        polynomials[indices[d]].value (p(d), tmp);
        v[d][0] = tmp[0];
        v[d][1] = tmp[1];
        v[d][2] = tmp[2];
      }
  }

  Tensor<2,dim> grad_grad;
  for (unsigned int d1=0; d1<dim; ++d1)
    for (unsigned int d2=0; d2<dim; ++d2)
      {
        grad_grad[d1][d2] = 1.;
        for (unsigned int x=0; x<dim; ++x)
          {
            unsigned int derivative=0;
            if (d1==x || d2==x)
              {
                if (d1==d2)
                  derivative=2;
                else
                  derivative=1;
              }
            grad_grad[d1][d2] *= v[x][derivative];
          }
      }

  return grad_grad;
}




template <int dim, typename POLY>
void
ABFScalarPolynomials<dim,POLY>::
compute (const Point<dim>            &p,
         std::vector<double>         &values,
         std::vector<Tensor<1,dim> > &grads,
         std::vector<Tensor<2,dim> > &grad_grads) const
{
  Assert (values.size()==n_tensor_pols    || values.size()==0,
          ExcDimensionMismatch2(values.size(), n_tensor_pols, 0));
  Assert (grads.size()==n_tensor_pols     || grads.size()==0,
          ExcDimensionMismatch2(grads.size(), n_tensor_pols, 0));
  Assert (grad_grads.size()==n_tensor_pols|| grad_grads.size()==0,
          ExcDimensionMismatch2(grad_grads.size(), n_tensor_pols, 0));

  const bool update_values     = (values.size() == n_tensor_pols),
             update_grads      = (grads.size()==n_tensor_pols),
             update_grad_grads = (grad_grads.size()==n_tensor_pols);

  // check how many
  // values/derivatives we have to
  // compute
  unsigned int n_values_and_derivatives = 0;
  if (update_values)
    n_values_and_derivatives = 1;
  if (update_grads)
    n_values_and_derivatives = 2;
  if (update_grad_grads)
    n_values_and_derivatives = 3;


  // compute the values (and derivatives, if
  // necessary) of all polynomials at this
  // evaluation point. to avoid many
  // reallocation, use one std::vector for
  // polynomial evaluation and store the
  // result as Tensor<1,3> (that has enough
  // fields for any evaluation of values and
  // derivatives)
  Table<2,Tensor<1,3> > v(dim, polynomials.size());
  {
    std::vector<double> tmp (n_values_and_derivatives);
    for (unsigned int d=0; d<dim; ++d)
      for (unsigned int i=0; i<polynomials.size(); ++i)
        {
          polynomials[i].value(p(d), tmp);
          for (unsigned int e=0; e<n_values_and_derivatives; ++e)
            v(d,i)[e] = tmp[e];
        };
  }

  for (unsigned int i=0; i<n_tensor_pols; ++i)
    {
      // first get the
      // one-dimensional indices of
      // this particular tensor
      // product polynomial
      unsigned int indices[dim];
      compute_index (i, indices);

      if (update_values)
        {
          values[i] = 1;
          for (unsigned int x=0; x<dim; ++x)
            values[i] *= v(x,indices[x])[0];
        }

      if (update_grads)
        for (unsigned int d=0; d<dim; ++d)
          {
            grads[i][d] = 1.;
            for (unsigned int x=0; x<dim; ++x)
              grads[i][d] *= v(x,indices[x])[d==x];
          }

      if (update_grad_grads)
        for (unsigned int d1=0; d1<dim; ++d1)
          for (unsigned int d2=0; d2<dim; ++d2)
            {
              grad_grads[i][d1][d2] = 1.;
              for (unsigned int x=0; x<dim; ++x)
                {
                  unsigned int derivative=0;
                  if (d1==x || d2==x)
                    {
                      if (d1==d2)
                        derivative=2;
                      else
                        derivative=1;
                    }
                  grad_grads[i][d1][d2]
                  *= v(x,indices[x])[derivative];
                }
            }
    }
}

template <int dim, typename POLY>
unsigned int
ABFScalarPolynomials<dim, POLY>::
get_n_tensor_pols (const std::vector<POLY>  &pols)
{
  Assert (pols.size()>=2, ExcMessage("1D Polynomial degree k<0"));
  Assert (dim>=1, ExcMessage("Dimension has to be >=1 "));\

  unsigned int y = 1;
  for (unsigned int d=0; d<dim; ++d)
    y *= pols.size() - 1;

  unsigned int x = 1;
  for (unsigned int d=0; d<dim-1; ++d)
    x *= pols.size() - 1;

  y += dim*x;

  return y;
}


/* ------------------- explicit instantiations -------------- */
template class ABFScalarPolynomials<1,Polynomials::Polynomial<double> >;
template class ABFScalarPolynomials<2,Polynomials::Polynomial<double> >;
template class ABFScalarPolynomials<3,Polynomials::Polynomial<double> >;

template class ABFScalarPolynomials<1,Polynomials::PiecewisePolynomial<double> >;
template class ABFScalarPolynomials<2,Polynomials::PiecewisePolynomial<double> >;
template class ABFScalarPolynomials<3,Polynomials::PiecewisePolynomial<double> >;

DEAL_II_NAMESPACE_CLOSE
