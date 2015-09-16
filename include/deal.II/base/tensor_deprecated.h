// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

#ifndef dealii__tensor_deprecated_h
#define dealii__tensor_deprecated_h

#include <deal.II/base/tensor.h>


/* --------- Deprecated non-member functions operating on tensors. ---------- */

DEAL_II_NAMESPACE_OPEN

/**
 * @name Deprecated Tensor operations
 */
//@{

/**
 * Single contraction for tensors: contract the last index of a tensor @p
 * src1 of rank @p rank_1 with the first index of a tensor @p src2 of rank
 * @p rank_2.
 *
 * @deprecated Use operator* instead. It denotes a single contraction.
 * @relates Tensor
 */
template <int rank_1, int rank_2, int dim, typename Number>
inline
void contract (Tensor<rank_1 + rank_2 - 2, dim, Number> &dest,
               const Tensor<rank_1 ,dim, Number>        &src1,
               const Tensor<rank_2 ,dim, Number>        &src2) DEAL_II_DEPRECATED;

/**
 * Contract a tensor of rank 1 with a tensor of rank 1 and return the
 * result.
 *
 * @deprecated Use operator* instead. It denotes a single contraction.
 * @relates Tensor
 */
template <int dim, typename Number, typename OtherNumber>
inline
typename ProductType<Number,OtherNumber>::type
contract (const Tensor<1,dim,Number> &src1,
          const Tensor<1,dim,OtherNumber> &src2) DEAL_II_DEPRECATED;

/**
 * The cross-product of 2 vectors in 3d.
 *
 * @deprecated Use the cross_product function that returns the value.
 * @relates Tensor
 */
template <int dim, typename Number>
inline
void
cross_product (Tensor<1,dim,Number>       &dst,
               const Tensor<1,dim,Number> &src1,
               const Tensor<1,dim,Number> &src2) DEAL_II_DEPRECATED;

/**
 * Form the outer product of two tensors.
 *
 * @deprecated Use the generic version that returns its result instead.
 * @relates Tensor
 */
template <int rank_1, int rank_2, int dim, typename Number>
void outer_product(Tensor<rank_1 + rank_2, dim, Number> &dst,
                   const Tensor<rank_1, dim, Number>    &src1,
                   const Tensor<rank_2, dim, Number>    &src2) DEAL_II_DEPRECATED;

/**
 * Multiply a Tensor<1,dim,Number> with a Number.
 *
 * @deprecated Use operator* instead.
 * @relates Tensor
 */
template <int dim, typename Number>
void outer_product (Tensor<1,dim,Number>       &dst,
                    const Number                src1,
                    const Tensor<1,dim,Number> &src2) DEAL_II_DEPRECATED;

/**
 * Multiply a Tensor<1,dim,Number> with a Number.
 *
 * @deprecated Use operator* instead.
 * @relates Tensor
 */
template <int dim, typename Number>
void outer_product (Tensor<1,dim,Number>       &dst,
                    const Tensor<1,dim,Number>  src1,
                    const Number                src2) DEAL_II_DEPRECATED;

/**
 * @deprecated Do not use this function, evaluate the value manually.
 * @relates Tensor
 */
template <int rank, typename Number>
inline
Number determinant (const Tensor<rank,1,Number> &t) DEAL_II_DEPRECATED;


/**
 * @deprecated Do not use this function, evaluate the value manually.
 * @relates Tensor
 */
template <typename Number>
inline
Number determinant (const Tensor<1,1,Number> &t) DEAL_II_DEPRECATED;


template <int rank_1, int rank_2, int dim, typename Number>
inline
void contract (Tensor<rank_1 + rank_2 - 2, dim, Number> &dest,
               const Tensor<rank_1 ,dim, Number>        &src1,
               const Tensor<rank_2 ,dim, Number>        &src2)
{
  TensorAccessors::internal::ReorderedIndexView<0, rank_2, const Tensor<rank_2, dim, Number> >
  reordered = TensorAccessors::reordered_index_view<0, rank_2>(src2);
  TensorAccessors::contract<1, rank_1, rank_2, dim>(dest, src1, reordered);
}

template <int dim, typename Number, typename OtherNumber>
inline
typename ProductType<Number,OtherNumber>::type
contract (const Tensor<1,dim,Number> &src1,
          const Tensor<1,dim,OtherNumber> &src2)
{
  typename ProductType<Number,OtherNumber>::type res
    = typename ProductType<Number,OtherNumber>::type();
  for (unsigned int i=0; i<dim; ++i)
    res += src1[i] * src2[i];

  return res;
}

template <int dim, typename Number>
inline
void
cross_product (Tensor<1,dim,Number>       &dst,
               const Tensor<1,dim,Number> &src1,
               const Tensor<1,dim,Number> &src2)
{
  dst = cross_product(src1, src2);
}

template <int rank_1, int rank_2, int dim, typename Number>
void outer_product(Tensor<rank_1 + rank_2, dim, Number> &dst,
                   const Tensor<rank_1, dim, Number>    &src1,
                   const Tensor<rank_2, dim, Number>    &src2)
{
  TensorAccessors::contract<0, rank_1, rank_2, dim>(dst, src1, src2);
}

template <int dim, typename Number>
void outer_product (Tensor<1,dim,Number>       &dst,
                    const Number                src1,
                    const Tensor<1,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    dst[i] = src1 * src2[i];
}

template <int dim, typename Number>
void outer_product (Tensor<1,dim,Number>       &dst,
                    const Tensor<1,dim,Number>  src1,
                    const Number         src2)
{
  for (unsigned int i=0; i<dim; ++i)
    dst[i] = src1[i] * src2;
}

template <int rank, typename Number>
inline
Number determinant (const Tensor<rank,1,Number> &t)
{
  return determinant(t[0]);
}

template <typename Number>
inline
Number determinant (const Tensor<1,1,Number> &t)
{
  return t[0];
}

//@}

DEAL_II_NAMESPACE_CLOSE

#endif
