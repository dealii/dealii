//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/tensor_product_polynomials.h>
#include <base/exceptions.h>
#include <base/table.h>



/* ------------------- TensorProductPolynomials -------------- */

template <>
void
TensorProductPolynomials<1>::
compute_index (const unsigned int i,
               unsigned int       (&indices)[1]) const
{
  Assert (i<polynomials.size(), ExcInternalError());
  indices[0] = index_map[i];
}



template <>
void
TensorProductPolynomials<2>::
compute_index (const unsigned int i,
               unsigned int       (&indices)[2]) const
{
  const unsigned int n_pols = polynomials.size();
  Assert (i<n_pols*n_pols, ExcInternalError());
  const unsigned int n=index_map[i];

  indices[0] = n % n_pols;
  indices[1] = n / n_pols;
}



template <>
void
TensorProductPolynomials<3>::
compute_index (const unsigned int i,
               unsigned int       (&indices)[3]) const
{
  const unsigned int n_pols = polynomials.size();
  Assert (i<n_pols*n_pols*n_pols, ExcInternalError());
  const unsigned int n=index_map[i];

  indices[0] = n % n_pols;
  indices[1] = (n/n_pols) % n_pols;
  indices[2] = n / (n_pols*n_pols);
}


template <int dim>
void
TensorProductPolynomials<dim>::output_indices(std::ostream &out) const
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



template <int dim>
void
TensorProductPolynomials<dim>::set_numbering(
  const std::vector<unsigned int> &renumber)
{
  Assert(renumber.size()==index_map.size(),
	 ExcDimensionMismatch(renumber.size(), index_map.size()));

  index_map=renumber;
  for (unsigned int i=0; i<index_map.size(); ++i)
    index_map_inverse[index_map[i]]=i;
}



template <int dim>
double
TensorProductPolynomials<dim>::compute_value (const unsigned int i,
                                              const Point<dim> &p) const
{
  unsigned int indices[dim];
  compute_index (i, indices);
  
  double value=1.;
  for (unsigned int d=0; d<dim; ++d)
    value *= polynomials[indices[d]].value(p(d));
  
  return value;
}

  
template <int dim>
Tensor<1,dim>
TensorProductPolynomials<dim>::compute_grad (const unsigned int i,
                                             const Point<dim> &p) const
{
  unsigned int indices[dim];
  compute_index (i, indices);

                                   // compute values and
                                   // uni-directional derivatives at
                                   // the given point in each
                                   // co-ordinate direction
  std::vector<std::vector<double> > v(dim, std::vector<double> (2));
  for (unsigned int d=0; d<dim; ++d)
    polynomials[indices[d]].value(p(d), v[d]);
  
  Tensor<1,dim> grad;
  for (unsigned int d=0; d<dim; ++d)
    {
      grad[d] = 1.;
      for (unsigned int x=0; x<dim; ++x)
        grad[d] *= v[x][d==x];
    }
  
  return grad;
}


template <int dim>
Tensor<2,dim>
TensorProductPolynomials<dim>::compute_grad_grad (const unsigned int i,
                                                  const Point<dim> &p) const
{
  unsigned int indices[dim];
  compute_index (i, indices);

  std::vector<std::vector<double> > v(dim, std::vector<double> (3));
  for (unsigned int d=0; d<dim; ++d)
    polynomials[indices[d]].value(p(d), v[d]);
  
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




template <int dim>
void
TensorProductPolynomials<dim>::
compute (const Point<dim>            &p,
         std::vector<double>         &values,
         std::vector<Tensor<1,dim> > &grads,
         std::vector<Tensor<2,dim> > &grad_grads) const
{
  Assert (values.size()==n_tensor_pols || values.size()==0,
          ExcDimensionMismatch2(values.size(), n_tensor_pols, 0));
  Assert (grads.size()==n_tensor_pols|| grads.size()==0,
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


                                   // compute the values (and
                                   // derivatives, if necessary) of
                                   // all polynomials at this
                                   // evaluation point
  Table<2,std::vector<double> > v(dim, polynomials.size());
  for (unsigned int d=0; d<dim; ++d)
    for (unsigned int i=0; i<polynomials.size(); ++i)
      {
        v(d,i).resize (n_values_and_derivatives, 0.);
        polynomials[i].value(p(d), v(d,i));
      };
  
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



template<int dim>
unsigned int
TensorProductPolynomials<dim>::n() const
{
  return n_tensor_pols;
}




/* ------------------- AnisotropicPolynomials -------------- */


template <int dim>
AnisotropicPolynomials<dim>::
AnisotropicPolynomials(const std::vector<std::vector<Polynomials::Polynomial<double> > > &pols)
		:
		polynomials (pols),
		n_tensor_pols(get_n_tensor_pols(pols))
{
  Assert (pols.size() == dim, ExcDimensionMismatch(pols.size(), dim));
  for (unsigned int d=0; d<dim; ++d)
    Assert (pols[d].size() > 0,
            ExcMessage ("The number of polynomials must be larger than zero "
                        "for all coordinate directions."));
}




template <>
void
AnisotropicPolynomials<1>::
compute_index (const unsigned int i,
               unsigned int       (&indices)[1]) const
{
  Assert (i<polynomials[0].size(), ExcInternalError());
  indices[0] = i;
}



template <>
void
AnisotropicPolynomials<2>::
compute_index (const unsigned int i,
               unsigned int       (&indices)[2]) const
{
  Assert (i < polynomials[0].size()*polynomials[1].size(),
          ExcInternalError());

  indices[0] = i % polynomials[0].size();
  indices[1] = i / polynomials[0].size();
}



template <>
void
AnisotropicPolynomials<3>::
compute_index (const unsigned int i,
               unsigned int       (&indices)[3]) const
{
  Assert (i < polynomials[0].size()*polynomials[1].size()*polynomials[2].size(),
          ExcInternalError());

  indices[0] = i % polynomials[0].size();
  indices[1] = (i/polynomials[0].size()) % polynomials[1].size();
  indices[2] = i / (polynomials[0].size()*polynomials[1].size());
}



template <int dim>
double
AnisotropicPolynomials<dim>::compute_value (const unsigned int i,
                                              const Point<dim> &p) const
{
  unsigned int indices[dim];
  compute_index (i, indices);
  
  double value=1.;
  for (unsigned int d=0; d<dim; ++d)
    value *= polynomials[d][indices[d]].value(p(d));
  
  return value;
}

  
template <int dim>
Tensor<1,dim>
AnisotropicPolynomials<dim>::compute_grad (const unsigned int i,
                                             const Point<dim> &p) const
{
  unsigned int indices[dim];
  compute_index (i, indices);

                                   // compute values and
                                   // uni-directional derivatives at
                                   // the given point in each
                                   // co-ordinate direction
  std::vector<std::vector<double> > v(dim, std::vector<double> (2));
  for (unsigned int d=0; d<dim; ++d)
    polynomials[d][indices[d]].value(p(d), v[d]);
  
  Tensor<1,dim> grad;
  for (unsigned int d=0; d<dim; ++d)
    {
      grad[d] = 1.;
      for (unsigned int x=0; x<dim; ++x)
        grad[d] *= v[x][d==x];
    }
  
  return grad;
}


template <int dim>
Tensor<2,dim>
AnisotropicPolynomials<dim>::compute_grad_grad (const unsigned int i,
                                                  const Point<dim> &p) const
{
  unsigned int indices[dim];
  compute_index (i, indices);

  std::vector<std::vector<double> > v(dim, std::vector<double> (3));
  for (unsigned int d=0; d<dim; ++d)
    polynomials[d][indices[d]].value(p(d), v[d]);
  
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




template <int dim>
void
AnisotropicPolynomials<dim>::
compute (const Point<dim>            &p,
         std::vector<double>         &values,
         std::vector<Tensor<1,dim> > &grads,
         std::vector<Tensor<2,dim> > &grad_grads) const
{
  Assert (values.size()==n_tensor_pols || values.size()==0,
          ExcDimensionMismatch2(values.size(), n_tensor_pols, 0));
  Assert (grads.size()==n_tensor_pols|| grads.size()==0,
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


                                   // compute the values (and
                                   // derivatives, if necessary) of
                                   // all polynomials at this
                                   // evaluation point
  std::vector<std::vector<std::vector<double> > > v(dim);
  for (unsigned int d=0; d<dim; ++d)
    {
      v[d].resize (polynomials[d].size());
      for (unsigned int i=0; i<polynomials[d].size(); ++i)
        {
          v[d][i].resize (n_values_and_derivatives, 0.);
          polynomials[d][i].value(p(d), v[d][i]);
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
            values[i] *= v[x][indices[x]][0];
        }
  
      if (update_grads)
        for (unsigned int d=0; d<dim; ++d)
          {
            grads[i][d] = 1.;            
            for (unsigned int x=0; x<dim; ++x)
              grads[i][d] *= v[x][indices[x]][d==x ? 1 : 0];
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
                    *= v[x][indices[x]][derivative];
                }
            }
    }
}



template<int dim>
unsigned int
AnisotropicPolynomials<dim>::n() const
{
  return n_tensor_pols;
}


template <int dim>
unsigned int
AnisotropicPolynomials<dim>::
get_n_tensor_pols (const std::vector<std::vector<Polynomials::Polynomial<double> > > &pols)
{
  unsigned int y = 1;
  for (unsigned int d=0; d<dim; ++d)
    y *= pols[d].size();
  return y;
}



/* ------------------- explicit instantiations -------------- */
template class TensorProductPolynomials<1>;
template class TensorProductPolynomials<2>;
template class TensorProductPolynomials<3>;

template class AnisotropicPolynomials<1>;
template class AnisotropicPolynomials<2>;
template class AnisotropicPolynomials<3>;
