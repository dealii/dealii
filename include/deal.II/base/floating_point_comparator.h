// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_base_floating_point_copmerator_h
#define dealii_base_floating_point_copmerator_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * A class that is used to compare floating point arrays (e.g. std::vector,
 * Tensor<1,dim>, etc.). The most common use case of this comparator
 * is for detecting arrays with the same content for the sake of
 * compression. The idea of this class is to consider two arrays as
 * equal if they are the same within a given tolerance. We use this
 * comparator class within a std::map<> of the given arrays. Note that this
 * comparison operator does not satisfy all the mathematical properties one
 * usually wants to have (consider e.g. the numbers a=0, b=0.1, c=0.2 with
 * tolerance 0.15; the operator gives a<c, but neither a<b? nor b<c? is
 * satisfied). This is not a problem in the use cases for this class, but be
 * careful when using it in other contexts.
 */
template <typename VectorizedArrayType>
struct FloatingPointComparator
{
  using Number = typename dealii::internal::VectorizedArrayTrait<
    VectorizedArrayType>::value_type;
  static constexpr std::size_t width =
    dealii::internal::VectorizedArrayTrait<VectorizedArrayType>::width;

  FloatingPointComparator(const Number scaling);

  /**
   * Compare two vectors of numbers (not necessarily of the same length),
   * where vectors of different lengths are first sorted by their length and
   * then by the entries.
   */
  bool
  operator()(const std::vector<Number> &v1,
             const std::vector<Number> &v2) const;

  /**
   * Compare two vectorized arrays (stored as tensors to avoid alignment
   * issues).
   */
  bool
  operator()(const Tensor<1, width, Number> &t1,
             const Tensor<1, width, Number> &t2) const;

  /**
   * Compare two rank-1 tensors of vectorized arrays (stored as tensors to
   * avoid alignment issues).
   */
  template <int dim>
  bool
  operator()(const Tensor<1, dim, Tensor<1, width, Number>> &t1,
             const Tensor<1, dim, Tensor<1, width, Number>> &t2) const;

  /**
   * Compare two rank-2 tensors of vectorized arrays (stored as tensors to
   * avoid alignment issues).
   */
  template <int dim>
  bool
  operator()(const Tensor<2, dim, Tensor<1, width, Number>> &t1,
             const Tensor<2, dim, Tensor<1, width, Number>> &t2) const;

  /**
   * Compare two arrays of tensors.
   */
  template <int dim>
  bool
  operator()(const std::array<Tensor<2, dim, Number>, dim + 1> &t1,
             const std::array<Tensor<2, dim, Number>, dim + 1> &t2) const;

  Number tolerance;
};


/* ------------------------------------------------------------------ */


template <typename VectorizedArrayType>
FloatingPointComparator<VectorizedArrayType>::FloatingPointComparator(
  const Number scaling)
  : tolerance(scaling * std::numeric_limits<double>::epsilon() * 1024.)
{}



template <typename VectorizedArrayType>
bool
FloatingPointComparator<VectorizedArrayType>::operator()(
  const std::vector<Number> &v1,
  const std::vector<Number> &v2) const
{
  const unsigned int s1 = v1.size(), s2 = v2.size();
  if (s1 < s2)
    return true;
  else if (s1 > s2)
    return false;
  else
    for (unsigned int i = 0; i < s1; ++i)
      if (v1[i] < v2[i] - tolerance)
        return true;
      else if (v1[i] > v2[i] + tolerance)
        return false;
  return false;
}



template <typename VectorizedArrayType>
bool
FloatingPointComparator<VectorizedArrayType>::operator()(
  const Tensor<1, width, Number> &t1,
  const Tensor<1, width, Number> &t2) const
{
  for (unsigned int k = 0; k < width; ++k)
    if (t1[k] < t2[k] - tolerance)
      return true;
    else if (t1[k] > t2[k] + tolerance)
      return false;
  return false;
}



template <typename VectorizedArrayType>
template <int dim>
bool
FloatingPointComparator<VectorizedArrayType>::operator()(
  const Tensor<1, dim, Tensor<1, width, Number>> &t1,
  const Tensor<1, dim, Tensor<1, width, Number>> &t2) const
{
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int k = 0; k < width; ++k)
      if (t1[d][k] < t2[d][k] - tolerance)
        return true;
      else if (t1[d][k] > t2[d][k] + tolerance)
        return false;
  return false;
}



template <typename VectorizedArrayType>
template <int dim>
bool
FloatingPointComparator<VectorizedArrayType>::operator()(
  const Tensor<2, dim, Tensor<1, width, Number>> &t1,
  const Tensor<2, dim, Tensor<1, width, Number>> &t2) const
{
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int e = 0; e < dim; ++e)
      for (unsigned int k = 0; k < width; ++k)
        if (t1[d][e][k] < t2[d][e][k] - tolerance)
          return true;
        else if (t1[d][e][k] > t2[d][e][k] + tolerance)
          return false;
  return false;
}



template <typename VectorizedArrayType>
template <int dim>
bool
FloatingPointComparator<VectorizedArrayType>::operator()(
  const std::array<Tensor<2, dim, Number>, dim + 1> &t1,
  const std::array<Tensor<2, dim, Number>, dim + 1> &t2) const
{
  for (unsigned int i = 0; i < t1.size(); ++i)
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int e = 0; e < dim; ++e)
        if (t1[i][d][e] < t2[i][d][e] - tolerance)
          return true;
        else if (t1[i][d][e] > t2[i][d][e] + tolerance)
          return false;
  return false;
}


DEAL_II_NAMESPACE_CLOSE

#endif
