// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_base_floating_point_copmerator_h
#define dealii_base_floating_point_copmerator_h

#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
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
    dealii::internal::VectorizedArrayTrait<Number>::width();

  /**
   * An enum to decode whether a particular comparison returns less, larger,
   * or equal. This is needed beyond the regular operator() because of
   * efficiency reasons in nested checks as supported by this class.
   */
  enum class ComparisonResult
  {
    less,
    equal,
    greater
  };

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
   * Compare two objects of the same type and return `true` in case the first
   * object is less than the second one. This function calls the compare()
   * functions below, so only types supported by compare() will work with this
   * function.
   */
  template <typename T>
  bool
  operator()(const T &object1, const T &object2) const;

  /**
   * Compare two vectors of numbers (not necessarily of the same length),
   * where vectors of different lengths are first sorted by their length and
   * then by the entries.
   */
  template <typename T>
  ComparisonResult
  compare(const std::vector<T> &v1, const std::vector<T> &v2) const;

  /**
   * Compare two arrays.
   */
  template <std::size_t dim, typename T>
  ComparisonResult
  compare(const std::array<T, dim> &t1, const std::array<T, dim> &t2) const;

  /**
   * Compare two tensors.
   */
  template <int rank, int dim, typename T>
  ComparisonResult
  compare(const Tensor<rank, dim, T> &t1, const Tensor<rank, dim, T> &t2) const;

  /**
   * Compare two derivative forms.
   */
  template <int rank, int dim, int spacedim, typename T>
  ComparisonResult
  compare(const DerivativeForm<rank, dim, spacedim, T> &t1,
          const DerivativeForm<rank, dim, spacedim, T> &t2) const;

  /**
   * Compare two tables.
   */
  template <typename T>
  ComparisonResult
  compare(const Table<2, T> &t1, const Table<2, T> &t2) const;

  /**
   * Compare two scalar numbers.
   */
  ComparisonResult
  compare(const ScalarNumber s1, const ScalarNumber s2) const;

  /**
   * Compare two VectorizedArray instances.
   */
  ComparisonResult
  compare(const VectorizedArray<ScalarNumber, width> &v1,
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
{
  Assert(mask.count() > 0, ExcMessage("No component selected"));
}



template <typename Number>
template <typename T>
bool
FloatingPointComparator<Number>::operator()(const T &object1,
                                            const T &object2) const
{
  return compare(object1, object2) == ComparisonResult::less;
}



template <typename Number>
template <typename T>
typename FloatingPointComparator<Number>::ComparisonResult
FloatingPointComparator<Number>::compare(const std::vector<T> &v1,
                                         const std::vector<T> &v2) const
{
  const unsigned int s1 = v1.size(), s2 = v2.size();
  if (s1 < s2)
    return ComparisonResult::less;
  else if (s1 > s2)
    return ComparisonResult::greater;
  else
    for (unsigned int i = 0; i < s1; ++i)
      {
        const ComparisonResult result = compare(v1[i], v2[i]);
        if (result != ComparisonResult::equal)
          return result;
      }
  return ComparisonResult::equal;
}



template <typename Number>
template <std::size_t dim, typename T>
typename FloatingPointComparator<Number>::ComparisonResult
FloatingPointComparator<Number>::compare(const std::array<T, dim> &t1,
                                         const std::array<T, dim> &t2) const
{
  for (unsigned int i = 0; i < t1.size(); ++i)
    {
      const ComparisonResult result = compare(t1[i], t2[i]);
      if (result != ComparisonResult::equal)
        return result;
    }
  return ComparisonResult::equal;
}



template <typename Number>
template <int rank, int dim, typename T>
typename FloatingPointComparator<Number>::ComparisonResult
FloatingPointComparator<Number>::compare(const Tensor<rank, dim, T> &t1,
                                         const Tensor<rank, dim, T> &t2) const
{
  for (unsigned int i = 0; i < dim; ++i)
    {
      const ComparisonResult result = compare(t1[i], t2[i]);
      if (result != ComparisonResult::equal)
        return result;
    }
  return ComparisonResult::equal;
}



template <typename Number>
template <int rank, int dim, int spacedim, typename T>
typename FloatingPointComparator<Number>::ComparisonResult
FloatingPointComparator<Number>::compare(
  const DerivativeForm<rank, dim, spacedim, T> &t1,
  const DerivativeForm<rank, dim, spacedim, T> &t2) const
{
  for (unsigned int i = 0; i < spacedim; ++i)
    {
      const ComparisonResult result = compare(t1[i], t2[i]);
      if (result != ComparisonResult::equal)
        return result;
    }
  return ComparisonResult::equal;
}



template <typename Number>
template <typename T>
typename FloatingPointComparator<Number>::ComparisonResult
FloatingPointComparator<Number>::compare(const Table<2, T> &t1,
                                         const Table<2, T> &t2) const
{
  AssertDimension(t1.size(0), t2.size(0));
  AssertDimension(t1.size(1), t2.size(1));

  for (unsigned int i = 0; i < t1.size(0); ++i)
    for (unsigned int j = 0; j < t1.size(1); ++j)
      {
        const ComparisonResult result = compare(t1[i][j], t2[i][j]);
        if (result != ComparisonResult::equal)
          return result;
      }
  return ComparisonResult::equal;
}



template <typename Number>
typename FloatingPointComparator<Number>::ComparisonResult
FloatingPointComparator<Number>::compare(const ScalarNumber s1,
                                         const ScalarNumber s2) const
{
  if (width == 1 || mask[0])
    {
      const ScalarNumber tolerance =
        use_absolute_tolerance ?
          this->tolerance :
          ((std::abs(s1) + std::abs(s2)) * this->tolerance);

      if (s1 < s2 - tolerance)
        return ComparisonResult::less;
      else if (s1 > s2 + tolerance)
        return ComparisonResult::greater;
    }

  return ComparisonResult::equal;
}

template <typename Number>
typename FloatingPointComparator<Number>::ComparisonResult
FloatingPointComparator<Number>::compare(
  const VectorizedArray<ScalarNumber, width> &v1,
  const VectorizedArray<ScalarNumber, width> &v2) const
{
  for (unsigned int i = 0; i < width; ++i)
    if (mask[i])
      {
        const ComparisonResult result = compare(v1[i], v2[i]);
        if (result != ComparisonResult::equal)
          return result;
      }

  return ComparisonResult::equal;
}


DEAL_II_NAMESPACE_CLOSE

#endif
