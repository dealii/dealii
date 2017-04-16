// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2015 by the deal.II authors
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


#include <deal.II/base/tensor_product_polynomials_const.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>

DEAL_II_NAMESPACE_OPEN



/* ------------------- TensorProductPolynomialsConst -------------- */


template <int dim>
double
TensorProductPolynomialsConst<dim>::compute_value (const unsigned int i,
                                                   const Point<dim> &p) const
{
  const unsigned int max_indices = this->n_tensor_pols;
  Assert (i<=max_indices, ExcInternalError());

  // treat the regular basis functions
  if (i<max_indices)
    return this->TensorProductPolynomials<dim>::compute_value(i,p);
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
  const unsigned int max_indices = this->n_tensor_pols;
  Assert (i<=max_indices, ExcInternalError());

  // treat the regular basis functions
  if (i<max_indices)
    return this->TensorProductPolynomials<dim>::compute_grad(i,p);
  else
    // this is for the constant function
    return Tensor<1,dim>();
}

template <int dim>
Tensor<2,dim>
TensorProductPolynomialsConst<dim>::compute_grad_grad (const unsigned int i,
                                                       const Point<dim> &p) const
{
  const unsigned int max_indices = this->n_tensor_pols;
  Assert (i<=max_indices, ExcInternalError());

  // treat the regular basis functions
  if (i<max_indices)
    return this->TensorProductPolynomials<dim>::compute_grad_grad(i,p);
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
         std::vector<Tensor<2,dim> > &grad_grads,
         std::vector<Tensor<3,dim> > &third_derivatives,
         std::vector<Tensor<4,dim> > &fourth_derivatives) const
{
  Assert (values.size()==this->n_tensor_pols+1 || values.size()==0,
          ExcDimensionMismatch2(values.size(), this->n_tensor_pols+1, 0));
  Assert (grads.size()==this->n_tensor_pols+1     || grads.size()==0,
          ExcDimensionMismatch2(grads.size(), this->n_tensor_pols+1, 0));
  Assert (grad_grads.size()==this->n_tensor_pols+1 || grad_grads.size()==0,
          ExcDimensionMismatch2(grad_grads.size(), this->n_tensor_pols+1, 0));
  Assert (third_derivatives.size()==this->n_tensor_pols+1 || third_derivatives.size()==0,
          ExcDimensionMismatch2(third_derivatives.size(), this->n_tensor_pols+1, 0));
  Assert (fourth_derivatives.size()==this->n_tensor_pols+1 || fourth_derivatives.size()==0,
          ExcDimensionMismatch2(fourth_derivatives.size(), this->n_tensor_pols+1, 0));

  // remove slot for const value, go into the base class compute method and
  // finally append the const value again
  bool do_values = false, do_grads = false, do_grad_grads = false;
  bool do_3rd_derivatives = false, do_4th_derivatives = false;
  if (values.empty() == false)
    {
      values.pop_back();
      do_values = true;
    }
  if (grads.empty() == false)
    {
      grads.pop_back();
      do_grads = true;
    }
  if (grad_grads.empty() == false)
    {
      grad_grads.pop_back();
      do_grad_grads = true;
    }
  if (third_derivatives.empty() == false)
    {
      third_derivatives.resize(this->n_tensor_pols);
      do_3rd_derivatives = true;
    }
  if (fourth_derivatives.empty() == false)
    {
      fourth_derivatives.resize(this->n_tensor_pols);
      do_4th_derivatives = true;
    }

  this->TensorProductPolynomials<dim>::compute(p, values, grads, grad_grads, third_derivatives, fourth_derivatives);

  //for dgq node: values =1, grads=0, grads_grads=0, third_derivatives=0, fourth_derivatives=0
  if (do_values)
    values.push_back(1.);
  if (do_grads)
    grads.emplace_back();
  if (do_grad_grads)
    grad_grads.emplace_back();
  if (do_3rd_derivatives)
    third_derivatives.emplace_back();
  if (do_4th_derivatives)
    fourth_derivatives.emplace_back();
}


/* ------------------- explicit instantiations -------------- */
template class TensorProductPolynomialsConst<1>;
template class TensorProductPolynomialsConst<2>;
template class TensorProductPolynomialsConst<3>;

DEAL_II_NAMESPACE_CLOSE
