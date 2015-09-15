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

//@}
/**
 * @name Deprecated Tensor operations
 */
//@{


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


#ifndef DOXYGEN


template <int dim, typename Number>
inline
void
cross_product (Tensor<1,dim,Number>       &dst,
               const Tensor<1,dim,Number> &src1,
               const Tensor<1,dim,Number> &src2)
{
  dst = cross_product(src1, src2);
}

#endif /* DOXYGEN */

//@}

DEAL_II_NAMESPACE_CLOSE

#endif
