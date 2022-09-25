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

#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include <bitset>
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
template <typename Number>
struct FloatingPointComparator
{
  using ScalarNumber =
    typename dealii::internal::VectorizedArrayTrait<Number>::value_type;
  static constexpr std::size_t width =
    dealii::internal::VectorizedArrayTrait<Number>::width;

  /**
   * Constructor.
   */
  FloatingPointComparator(
    const ScalarNumber        tolerance,
    const bool                use_absolute_tolerance = true,
    const std::bitset<width> &mask = std::bitset<width>().flip());

  FloatingPointComparator(const FloatingPointComparator &rhs) = default;

  FloatingPointComparator(FloatingPointComparator &&rhs) noexcept = default;

  /**
   * Compare two vectors of numbers (not necessarily of the same length),
   * where vectors of different lengths are first sorted by their length and
   * then by the entries.
   */
  template <typename T>
  bool
  operator()(const std::vector<T> &v1, const std::vector<T> &v2) const;

  /**
   * Compare two arrays.
   */
  template <std::size_t dim, typename T>
  bool
  operator()(const std::array<T, dim> &t1, const std::array<T, dim> &t2) const;

  /**
   * Compare two tensors.
   */
  template <int rank, int dim, typename T>
  bool
  operator()(const Tensor<rank, dim, T> &t1,
             const Tensor<rank, dim, T> &t2) const;

  /**
   * Compare two tables.
   */
  template <typename T>
  bool
  operator()(const Table<2, T> &t1, const Table<2, T> &t2) const;

  /**
   * Compare two scalar numbers.
   */
  bool
  operator()(const ScalarNumber &s1, const ScalarNumber &s2) const;

  /**
   * Compare two  VectorizedArray instances.
   */
  bool
  operator()(const VectorizedArray<ScalarNumber, width> &v1,
             const VectorizedArray<ScalarNumber, width> &v2) const;

private:
  const ScalarNumber       tolerance;
  const bool               use_absolute_tolerance;
  const std::bitset<width> mask;
};


/* ------------------------------------------------------------------ */


template <typename Number>
FloatingPointComparator<Number>::FloatingPointComparator(
  const ScalarNumber        tolerance,
  const bool                use_absolute_tolerance,
  const std::bitset<width> &mask)
  : tolerance(tolerance)
  , use_absolute_tolerance(use_absolute_tolerance)
  , mask(mask)
{}



template <typename Number>
template <typename T>
bool
FloatingPointComparator<Number>::operator()(const std::vector<T> &v1,
                                            const std::vector<T> &v2) const
{
  const unsigned int s1 = v1.size(), s2 = v2.size();
  if (s1 < s2)
    return true;
  else if (s1 > s2)
    return false;
  else
    for (unsigned int i = 0; i < s1; ++i)
      if (this->operator()(v1[i], v2[i]))
        return true;
      else if (this->operator()(v2[i], v1[i]))
        return false;
  return false;
}



template <typename Number>
template <std::size_t dim, typename T>
bool
FloatingPointComparator<Number>::operator()(const std::array<T, dim> &t1,
                                            const std::array<T, dim> &t2) const
{
  for (unsigned int i = 0; i < t1.size(); ++i)
    if (this->operator()(t1[i], t2[i]))
      return true;
    else if (this->operator()(t2[i], t1[i]))
      return false;
  return false;
}



template <typename Number>
template <int rank, int dim, typename T>
bool
FloatingPointComparator<Number>::operator()(
  const Tensor<rank, dim, T> &t1,
  const Tensor<rank, dim, T> &t2) const
{
  for (unsigned int k = 0; k < dim; ++k)
    if (this->operator()(t1[k], t2[k]))
      return true;
    else if (this->operator()(t2[k], t1[k]))
      return false;
  return false;
}



template <typename Number>
template <typename T>
bool
FloatingPointComparator<Number>::operator()(const Table<2, T> &t1,
                                            const Table<2, T> &t2) const
{
  AssertDimension(t1.size(0), t2.size(0));
  AssertDimension(t1.size(1), t2.size(1));

  for (unsigned int i = 0; i < t1.size(0); ++i)
    for (unsigned int j = 0; j < t1.size(1); ++j)
      if (this->operator()(t1[i][j], t2[i][j]))
        return true;
      else if (this->operator()(t2[i][j], t1[i][j]))
        return false;
  return false;
}

template <typename Number>
bool
FloatingPointComparator<Number>::operator()(const ScalarNumber &s1,
                                            const ScalarNumber &s2) const
{
  if (mask[0])
    {
      const ScalarNumber tolerance = use_absolute_tolerance ?
                                       this->tolerance :
                                       (std::abs(s1 + s2) * this->tolerance);

      if ((s1 < s2 - tolerance))
        return true;
      else
        return false;
    }

  return false;
}

template <typename Number>
bool
FloatingPointComparator<Number>::operator()(
  const VectorizedArray<ScalarNumber, width> &v1,
  const VectorizedArray<ScalarNumber, width> &v2) const
{
  for (unsigned int v = 0; v < width; ++v)
    if (mask[v])
      {
        const ScalarNumber tolerance =
          use_absolute_tolerance ? this->tolerance :
                                   (std::abs(v1[v] + v2[v]) * this->tolerance);

        if (v1[v] < v2[v] - tolerance)
          return true;
        if (v1[v] > v2[v] + tolerance)
          return false;
      }

  return false;
}


DEAL_II_NAMESPACE_CLOSE

#endif
