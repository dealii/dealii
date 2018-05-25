// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#ifndef dealii_tensor_deprecated_h
#define dealii_tensor_deprecated_h

#include <deal.II/base/tensor.h>


/* --------- Deprecated non-member functions operating on tensors. ---------- */

DEAL_II_NAMESPACE_OPEN

/**
 * @name Deprecated Tensor operations
 */
//@{

/**
 * Exception.
 *
 * @deprecated
 */
DeclException1(
  ExcInvalidTensorContractionIndex,
  int,
  << "You have requested contraction of tensors over index " << arg1
  << ", but this is not possible for tensors of the current type.");


/**
 * Double contract two tensors of rank 2, thus computing the Frobenius inner
 * product <tt>sum<sub>i,j</sub> src1[i][j]*src2[i][j]</tt>.
 *
 * @deprecated Use the double_contract() function that takes indices as
 * template arguments and returns its result instead.
 * @relatesalso Tensor
 */
template <int dim, typename Number>
DEAL_II_DEPRECATED inline Number
double_contract(const Tensor<2, dim, Number> &src1,
                const Tensor<2, dim, Number> &src2);


/**
 * Contract the last two indices of <tt>src1</tt> with the two indices
 * <tt>src2</tt>, creating a rank-2 tensor. This is the matrix-vector product
 * analog operation between tensors of rank 4 and rank 2.
 *
 * @deprecated Use the double_contract() function that takes indices as
 * template arguments and returns its result instead.
 * @relatesalso Tensor
 */
template <int dim, typename Number>
DEAL_II_DEPRECATED inline void
  double_contract(Tensor<2, dim, Number> &      dest,
                  const Tensor<4, dim, Number> &src1,
                  const Tensor<2, dim, Number> &src2);

/**
 * Contract a tensor of rank 2 with a tensor of rank 2. The contraction is
 * performed over index <tt>index1</tt> of the first tensor, and
 * <tt>index2</tt> of the second tensor. Note that the number of the index is
 * counted from 1 on, not from zero as usual.
 *
 * @deprecated Use the contract() function that takes indices as template
 * arguments and returns its result instead.
 * @relatesalso Tensor
 */
template <int dim, typename Number>
DEAL_II_DEPRECATED inline void contract(Tensor<2, dim, Number> &      dest,
                                        const Tensor<2, dim, Number> &src1,
                                        const unsigned int            index1,
                                        const Tensor<2, dim, Number> &src2,
                                        const unsigned int            index3);

/**
 * Contract a tensor of rank 3 with a tensor of rank 1. The contraction is
 * performed over index <tt>index1</tt> of the first tensor. Note that the
 * number of the index is counted from 1 on, not from zero as usual.
 *
 * @deprecated Use the contract() function that takes indices as template
 * arguments and returns its result instead.
 * @relatesalso Tensor
 */
template <int dim, typename Number>
DEAL_II_DEPRECATED inline void contract(Tensor<2, dim, Number> &      dest,
                                        const Tensor<3, dim, Number> &src1,
                                        const unsigned int            index1,
                                        const Tensor<1, dim, Number> &src2);

/**
 * Contract a tensor of rank 3 with a tensor of rank 2. The contraction is
 * performed over index <tt>index1</tt> of the first tensor, and
 * <tt>index2</tt> of the second tensor. Note that the number of the index is
 * counted from 1 on, not from zero as usual.
 *
 * @deprecated Use the contract() function that takes indices as template
 * arguments and returns its result instead.
 * @relatesalso Tensor
 */
template <int dim, typename Number>
DEAL_II_DEPRECATED inline void contract(Tensor<3, dim, Number> &      dest,
                                        const Tensor<3, dim, Number> &src1,
                                        const unsigned int            index1,
                                        const Tensor<2, dim, Number> &src2,
                                        const unsigned int            index2);

/**
 * Single contraction for tensors: contract the last index of a tensor @p src1
 * of rank @p rank_1 with the first index of a tensor @p src2 of rank @p
 * rank_2.
 *
 * @deprecated Use operator* instead. It denotes a single contraction.
 * @relatesalso Tensor
 */
template <int rank_1, int rank_2, int dim, typename Number>
DEAL_II_DEPRECATED inline void
  contract(Tensor<rank_1 + rank_2 - 2, dim, Number> &dest,
           const Tensor<rank_1, dim, Number> &       src1,
           const Tensor<rank_2, dim, Number> &       src2);

/**
 * Contract a tensor of rank 1 with a tensor of rank 1 and return the result.
 *
 * @deprecated Use operator* instead. It denotes a single contraction.
 * @relatesalso Tensor
 */
template <int dim, typename Number, typename OtherNumber>
DEAL_II_DEPRECATED inline typename ProductType<Number, OtherNumber>::type
contract(const Tensor<1, dim, Number> &     src1,
         const Tensor<1, dim, OtherNumber> &src2);


/**
 * The cross product of one vector in 2d. This is just a rotation by 90
 * degrees.
 *
 * @deprecated Use the function cross_product_2d() that returns the value.
 * @relatesalso Tensor
 */
template <int dim, typename Number>
DEAL_II_DEPRECATED inline void cross_product(Tensor<1, dim, Number> &      dst,
                                             const Tensor<1, dim, Number> &src);

/**
 * The cross product of 2 vectors in 3d.
 *
 * @deprecated Use the function cross_product_3d() that returns the value.
 * @relatesalso Tensor
 */
template <int dim, typename Number>
DEAL_II_DEPRECATED inline void
  cross_product(Tensor<1, dim, Number> &      dst,
                const Tensor<1, dim, Number> &src1,
                const Tensor<1, dim, Number> &src2);

/**
 * Form the outer product of two tensors.
 *
 * @deprecated Use the generic version that returns its result instead.
 * @relatesalso Tensor
 */
template <int rank_1, int rank_2, int dim, typename Number>
DEAL_II_DEPRECATED inline void
outer_product(Tensor<rank_1 + rank_2, dim, Number> &dst,
              const Tensor<rank_1, dim, Number> &   src1,
              const Tensor<rank_2, dim, Number> &   src2);

/**
 * Multiply a Tensor<1,dim,Number> with a Number.
 *
 * @deprecated Use operator* instead.
 * @relatesalso Tensor
 */
template <int dim, typename Number>
DEAL_II_DEPRECATED inline void
  outer_product(Tensor<1, dim, Number> &      dst,
                const Number                  src1,
                const Tensor<1, dim, Number> &src2);

/**
 * Multiply a Tensor<1,dim,Number> with a Number.
 *
 * @deprecated Use operator* instead.
 * @relatesalso Tensor
 */
template <int dim, typename Number>
DEAL_II_DEPRECATED inline void outer_product(Tensor<1, dim, Number> &     dst,
                                             const Tensor<1, dim, Number> src1,
                                             const Number                 src2);

/**
 * @deprecated Do not use this function, evaluate the value manually.
 * @relatesalso Tensor
 */
template <int rank, typename Number>
DEAL_II_DEPRECATED inline Number
determinant(const Tensor<rank, 1, Number> &t);


/**
 * @deprecated Do not use this function, evaluate the value manually.
 * @relatesalso Tensor
 */
template <typename Number>
DEAL_II_DEPRECATED inline Number
determinant(const Tensor<1, 1, Number> &t);

//@}

/* ----------------------------- Definitions: ------------------------------- */

template <int dim, typename Number>
inline Number
double_contract(const Tensor<2, dim, Number> &src1,
                const Tensor<2, dim, Number> &src2)
{
  Number res = internal::NumberType<Number>::value(0.0);
  for (unsigned int i = 0; i < dim; ++i)
    res += src1[i] * src2[i];

  return res;
}

template <int dim, typename Number>
inline void double_contract(Tensor<2, dim, Number> &      dest,
                            const Tensor<4, dim, Number> &src1,
                            const Tensor<2, dim, Number> &src2)
{
  dest.clear();
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          dest[i][j] += src1[i][j][k][l] * src2[k][l];
}

template <int dim, typename Number>
inline void contract(Tensor<2, dim, Number> &      dest,
                     const Tensor<2, dim, Number> &src1,
                     const unsigned int            index1,
                     const Tensor<2, dim, Number> &src2,
                     const unsigned int            index2)
{
  dest.clear();

  switch (index1)
    {
      case 1:
        switch (index2)
          {
            case 1:
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  for (unsigned int k = 0; k < dim; ++k)
                    dest[i][j] += src1[k][i] * src2[k][j];
              break;
            case 2:
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  for (unsigned int k = 0; k < dim; ++k)
                    dest[i][j] += src1[k][i] * src2[j][k];
              break;

            default:
              Assert(false, (ExcInvalidTensorContractionIndex(index2)));
          };
        break;
      case 2:
        switch (index2)
          {
            case 1:
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  for (unsigned int k = 0; k < dim; ++k)
                    dest[i][j] += src1[i][k] * src2[k][j];
              break;
            case 2:
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  for (unsigned int k = 0; k < dim; ++k)
                    dest[i][j] += src1[i][k] * src2[j][k];
              break;

            default:
              Assert(false, (ExcInvalidTensorContractionIndex(index2)));
          };
        break;

      default:
        Assert(false, (ExcInvalidTensorContractionIndex(index1)));
    };
}

template <int dim, typename Number>
inline void contract(Tensor<2, dim, Number> &      dest,
                     const Tensor<3, dim, Number> &src1,
                     const unsigned int            index1,
                     const Tensor<1, dim, Number> &src2)
{
  dest.clear();

  switch (index1)
    {
      case 1:
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            for (unsigned int k = 0; k < dim; ++k)
              dest[i][j] += src1[k][i][j] * src2[k];
        break;

      case 2:
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            for (unsigned int k = 0; k < dim; ++k)
              dest[i][j] += src1[i][k][j] * src2[k];
        break;

      case 3:
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            for (unsigned int k = 0; k < dim; ++k)
              dest[i][j] += src1[i][j][k] * src2[k];
        break;

      default:
        Assert(false, (ExcInvalidTensorContractionIndex(index1)));
    };
}

template <int dim, typename Number>
inline void contract(Tensor<3, dim, Number> &      dest,
                     const Tensor<3, dim, Number> &src1,
                     const unsigned int            index1,
                     const Tensor<2, dim, Number> &src2,
                     const unsigned int            index2)
{
  dest.clear();

  switch (index1)
    {
      case 1:
        switch (index2)
          {
            case 1:
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  for (unsigned int k = 0; k < dim; ++k)
                    for (unsigned int l = 0; l < dim; ++l)
                      dest[i][j][k] += src1[l][i][j] * src2[l][k];
              break;
            case 2:
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  for (unsigned int k = 0; k < dim; ++k)
                    for (unsigned int l = 0; l < dim; ++l)
                      dest[i][j][k] += src1[l][i][j] * src2[k][l];
              break;
            default:
              Assert(false, (ExcInvalidTensorContractionIndex(index2)));
          }

        break;
      case 2:
        switch (index2)
          {
            case 1:
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  for (unsigned int k = 0; k < dim; ++k)
                    for (unsigned int l = 0; l < dim; ++l)
                      dest[i][j][k] += src1[i][l][j] * src2[l][k];
              break;
            case 2:
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  for (unsigned int k = 0; k < dim; ++k)
                    for (unsigned int l = 0; l < dim; ++l)
                      dest[i][j][k] += src1[i][l][j] * src2[k][l];
              break;
            default:
              Assert(false, (ExcInvalidTensorContractionIndex(index2)));
          }

        break;
      case 3:
        switch (index2)
          {
            case 1:
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  for (unsigned int k = 0; k < dim; ++k)
                    for (unsigned int l = 0; l < dim; ++l)
                      dest[i][j][k] += src1[i][j][l] * src2[l][k];
              break;
            case 2:
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  for (unsigned int k = 0; k < dim; ++k)
                    for (unsigned int l = 0; l < dim; ++l)
                      dest[i][j][k] += src1[i][j][l] * src2[k][l];
              break;
            default:
              Assert(false, (ExcInvalidTensorContractionIndex(index2)));
          }

        break;
      default:
        Assert(false, (ExcInvalidTensorContractionIndex(index1)));
    }
}

template <int rank_1, int rank_2, int dim, typename Number>
inline void contract(Tensor<rank_1 + rank_2 - 2, dim, Number> &dest,
                     const Tensor<rank_1, dim, Number> &       src1,
                     const Tensor<rank_2, dim, Number> &       src2)
{
  TensorAccessors::internal::
    ReorderedIndexView<0, rank_2, const Tensor<rank_2, dim, Number>>
      reordered = TensorAccessors::reordered_index_view<0, rank_2>(src2);
  TensorAccessors::contract<1, rank_1, rank_2, dim>(dest, src1, reordered);
}

template <int dim, typename Number, typename OtherNumber>
inline typename ProductType<Number, OtherNumber>::type
contract(const Tensor<1, dim, Number> &     src1,
         const Tensor<1, dim, OtherNumber> &src2)
{
  typename ProductType<Number, OtherNumber>::type res =
    typename ProductType<Number, OtherNumber>::type();
  for (unsigned int i = 0; i < dim; ++i)
    res += src1[i] * src2[i];

  return res;
}

template <int dim, typename Number>
inline void cross_product(Tensor<1, dim, Number> &      dst,
                          const Tensor<1, dim, Number> &src)
{
  dst = cross_product_2d(src);
}

template <int dim, typename Number>
inline void cross_product(Tensor<1, dim, Number> &      dst,
                          const Tensor<1, dim, Number> &src1,
                          const Tensor<1, dim, Number> &src2)
{
  dst = cross_product_3d(src1, src2);
}

template <int rank_1, int rank_2, int dim, typename Number>
inline void
outer_product(Tensor<rank_1 + rank_2, dim, Number> &dst,
              const Tensor<rank_1, dim, Number> &   src1,
              const Tensor<rank_2, dim, Number> &   src2)
{
  TensorAccessors::contract<0, rank_1, rank_2, dim>(dst, src1, src2);
}

template <int dim, typename Number>
inline void outer_product(Tensor<1, dim, Number> &      dst,
                          const Number                  src1,
                          const Tensor<1, dim, Number> &src2)
{
  for (unsigned int i = 0; i < dim; ++i)
    dst[i] = src1 * src2[i];
}

template <int dim, typename Number>
inline void outer_product(Tensor<1, dim, Number> &     dst,
                          const Tensor<1, dim, Number> src1,
                          const Number                 src2)
{
  for (unsigned int i = 0; i < dim; ++i)
    dst[i] = src1[i] * src2;
}

template <int rank, typename Number>
inline Number
determinant(const Tensor<rank, 1, Number> &t)
{
  return determinant(t[0]);
}

template <typename Number>
inline Number
determinant(const Tensor<1, 1, Number> &t)
{
  return t[0];
}

DEAL_II_NAMESPACE_CLOSE

#endif
