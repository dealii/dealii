// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_table_indices_h
#define dealii_table_indices_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>

#include <algorithm>
#include <iterator>
#include <ostream>


DEAL_II_NAMESPACE_OPEN


/**
 * A class representing a fixed size array of indices.
 *
 * It is used in tensorial objects like the TableBase and SymmetricTensor
 * classes to represent a nested choice of indices.
 *
 * @tparam N The number of indices stored in each object.
 *
 * @ingroup data
 */
template <int N>
class TableIndices
{
public:
  static_assert(N > 0,
                "TableIndices objects need to represent at least one index.");


  /**
   * Default constructor. This constructor sets all indices to zero.
   */
  constexpr TableIndices() = default;

  /**
   * Constructor. Initializes the indices stored by this object by the given
   * arguments @p indices
   *
   * This constructor will result in a compiler error if
   * the number of arguments given is different from the number of the
   * indices this class stores (i.e., the template argument `N` of
   * this class), or if any of the arguments is not of some integer type.
   */
  template <typename... T>
  constexpr TableIndices(const T... indices);

  /**
   * Read-only access the value of the <tt>i</tt>th index.
   */
  constexpr std::size_t
  operator[](const unsigned int i) const;

  /**
   * Write access the value of the <tt>i</tt>th index.
   */
  constexpr std::size_t &
  operator[](const unsigned int i);

  /**
   * Compare two index fields for equality.
   */
  constexpr bool
  operator==(const TableIndices<N> &other) const;

  /**
   * Compare two index fields for inequality.
   */
  constexpr bool
  operator!=(const TableIndices<N> &other) const;

  /**
   * Sort the indices in ascending order. While this operation is not very
   * useful for Table objects, it is used for the SymmetricTensor class.
   */
  DEAL_II_HOST
  DEAL_II_CONSTEXPR void
  sort();

  /**
   * Write or read the data of this object to or from a stream for the purpose
   * of serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

protected:
  /**
   * Store the indices in an array.
   */
  std::size_t indices[N]{};
};



/* --------------------- Template and inline functions ---------------- */

template <int N>
template <typename... T>
constexpr TableIndices<N>::TableIndices(const T... args)
  : indices{static_cast<std::size_t>(args)...}
{
  static_assert(
    internal::TemplateConstraints::all_true<std::is_integral_v<T>...>::value,
    "Not all of the parameters have integral type!");
  static_assert(sizeof...(T) == N, "Wrong number of constructor arguments!");
}


template <int N>
constexpr inline std::size_t
TableIndices<N>::operator[](const unsigned int i) const
{
  AssertIndexRange(i, N);
  return indices[i];
}


template <int N>
constexpr inline std::size_t &
TableIndices<N>::operator[](const unsigned int i)
{
  AssertIndexRange(i, N);
  return indices[i];
}


template <int N>
constexpr bool
TableIndices<N>::operator==(const TableIndices<N> &other) const
{
  return std::equal(std::begin(indices),
                    std::end(indices),
                    std::begin(other.indices));
}


template <int N>
constexpr bool
TableIndices<N>::operator!=(const TableIndices<N> &other) const
{
  return !(*this == other);
}


template <int N>
DEAL_II_HOST DEAL_II_CONSTEXPR inline void
TableIndices<N>::sort()
{
  std::sort(std::begin(indices), std::end(indices));
}


template <int N>
template <class Archive>
inline void
TableIndices<N>::serialize(Archive &ar, const unsigned int)
{
  ar &indices;
}


/**
 * Output operator for TableIndices objects; reports them in a list like this:
 * <code>[i1,i2,...]</code>.
 *
 * @relatesalso TableIndices
 */
template <int N>
std::ostream &
operator<<(std::ostream &out, const TableIndices<N> &indices)
{
  out << '[';
  for (unsigned int i = 0; i < N; ++i)
    {
      out << indices[i];
      if (i + 1 != N)
        out << ',';
    }
  out << ']';

  return out;
}


DEAL_II_NAMESPACE_CLOSE

#endif
