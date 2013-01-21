//---------------------------------------------------------------------------
//    $Id: trilinos_sparse_matrix.cc 27985 2013-01-08 21:36:23Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <deal.II/base/tensor_product_polynomials_const.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>

DEAL_II_NAMESPACE_OPEN



/* ------------------- TensorProductPolynomialsConst -------------- */

template <>
inline
void
TensorProductPolynomialsConst<1>::
compute_index (const unsigned int i,
               unsigned int       (&indices)[1]) const
{
  Assert (i<polynomials.size(), ExcInternalError());
  indices[0] = index_map[i];
}

template <>
inline
void
TensorProductPolynomialsConst<2>::
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
inline
void
TensorProductPolynomialsConst<3>::
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
TensorProductPolynomialsConst<dim>::output_indices(std::ostream &out) const
{
  unsigned int ix[dim];
  //ouput without constant locally constant function
  for (unsigned int i=0; i<n_tensor_pols-1; ++i)
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
TensorProductPolynomialsConst<dim>::set_numbering(
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
TensorProductPolynomialsConst<dim>::compute_value (const unsigned int i,
                                                   const Point<dim> &p) const
{
  const unsigned int max_indices = Utilities::fixed_power<dim>(polynomials.size());

  Assert (i<=max_indices, ExcInternalError());

  // treat the regular basis functions
  if (i<max_indices)
    {
      unsigned int indices[dim];
      compute_index (i, indices);

      double value=1.;
      for (unsigned int d=0; d<dim; ++d)
        value *= polynomials[indices[d]].value(p(d));

      return value;
    }
  else
    // this is for the constant function
    return 1.;
}

template <>
double
TensorProductPolynomialsConst<0>::compute_value (const unsigned int ,
                                                 const Point<0> &) const
{
  Assert (false, ExcNotImplemented());
  return 0.;
}


template <int dim>
Tensor<1,dim>
TensorProductPolynomialsConst<dim>::compute_grad (const unsigned int i,
                                                  const Point<dim> &p) const
{
  const unsigned int max_indices = Utilities::fixed_power<dim>(polynomials.size());

  Assert (i<=max_indices, ExcInternalError());

  // treat the regular basis functions
  if (i<max_indices)
    {
      Tensor<1,dim> grad;

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

      for (unsigned int d=0; d<dim; ++d)
        {
          grad[d] = 1.;
          for (unsigned int x=0; x<dim; ++x)
            grad[d] *= v[x][d==x];
        }

      return grad; //return 0 for q0 node
    }
  else
    // this is for the constant function
    return Tensor<1,dim>();
}



template <int dim>
Tensor<2,dim>
TensorProductPolynomialsConst<dim>::compute_grad_grad (const unsigned int i,
                                                       const Point<dim> &p) const
{
  const unsigned int max_indices = Utilities::fixed_power<dim>(polynomials.size());

  Assert (i<=max_indices, ExcInternalError());

  // treat the regular basis functions
  if (i<max_indices)
    {
      Tensor<2,dim> grad_grad;

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
      return grad_grad; //return 0 for q0 node
    }
  else
    // this is for the constant function
    return Tensor<2,dim>();
}

template <int dim>
void
TensorProductPolynomialsConst<dim>::
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
  Table<2,Tensor<1,3> > v(dim, polynomials.size()+1);
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

  for (unsigned int i=0; i<n_tensor_pols-1; ++i)
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
  //for dgq node: values =1, grads=0, grads_grads=0
  if (update_values)
    values[n_tensor_pols-1]=1;
}



template <int dim>
unsigned int
TensorProductPolynomialsConst<dim>::n() const
{
  return n_tensor_pols;
}



template <>
unsigned int
TensorProductPolynomialsConst<0>::n() const
{
  return numbers::invalid_unsigned_int;
}

/* ------------------- explicit instantiations -------------- */
template class TensorProductPolynomialsConst<1>;
template class TensorProductPolynomialsConst<2>;
template class TensorProductPolynomialsConst<3>;

DEAL_II_NAMESPACE_CLOSE
