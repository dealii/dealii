// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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

#ifndef dealii_table_indices_h
#define dealii_table_indices_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

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
 * @author Wolfgang Bangerth, Matthias Maier, 2002, 2015
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
   * Constructor. This is the appropriate constructor for an
   * object of type TableIndices<1> and initializes the single
   * index with @p index0.
   *
   * This constructor will result in a compiler error if
   * the template argument @p N is different from one.
   */
  constexpr explicit TableIndices(const std::size_t index0);

  /**
   * Constructor. This is the appropriate constructor for an
   * object of type TableIndices<2> and initializes the
   * indices stored by this object by the given arguments.
   *
   * This constructor will result in a compiler error if
   * the template argument @p N is different from two.
   */
  constexpr TableIndices(const std::size_t index0, const std::size_t index1);

  /**
   * Constructor. This is the appropriate constructor for an
   * object of type TableIndices<3> and initializes the
   * indices stored by this object by the given arguments.
   *
   * This constructor will result in a compiler error if
   * the template argument @p N is different from three.
   */
  constexpr TableIndices(const std::size_t index0,
                         const std::size_t index1,
                         const std::size_t index2);

  /**
   * Constructor. This is the appropriate constructor for an
   * object of type TableIndices<4> and initializes the
   * indices stored by this object by the given arguments.
   *
   * This constructor will result in a compiler error if
   * the template argument @p N is different from four.
   */
  constexpr TableIndices(const std::size_t index0,
                         const std::size_t index1,
                         const std::size_t index2,
                         const std::size_t index3);

  /**
   * Constructor. This is the appropriate constructor for an
   * object of type TableIndices<5> and initializes the
   * indices stored by this object by the given arguments.
   *
   * This constructor will result in a compiler error if
   * the template argument @p N is different from five.
   */
  constexpr TableIndices(const std::size_t index0,
                         const std::size_t index1,
                         const std::size_t index2,
                         const std::size_t index3,
                         const std::size_t index4);

  /**
   * Convenience constructor that takes up to 9 arguments. It can be used to
   * populate a TableIndices object upon creation, either completely, or
   * partially.
   *
   * Index entries that are not set by these arguments (either because they
   * are omitted, or because $N > 9$) are set to
   * numbers::invalid_unsigned_int.
   *
   * Note that only the first <tt>N</tt> arguments are actually used.
   *
   * @deprecated Use the constructor with the appropriate number of arguments
   *   to initialize the @p N indices instead.
   */
  DEAL_II_DEPRECATED
  TableIndices(const std::size_t index0,
               const std::size_t index1,
               const std::size_t index2,
               const std::size_t index3,
               const std::size_t index4,
               const std::size_t index5,
               const std::size_t index6 = numbers::invalid_unsigned_int,
               const std::size_t index7 = numbers::invalid_unsigned_int,
               const std::size_t index8 = numbers::invalid_unsigned_int);

  /**
   * Read-only access the value of the <tt>i</tt>th index.
   */
  DEAL_II_CONSTEXPR std::size_t operator[](const unsigned int i) const;

  /**
   * Write access the value of the <tt>i</tt>th index.
   */
  DEAL_II_CONSTEXPR std::size_t &operator[](const unsigned int i);

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
  DEAL_II_CONSTEXPR void
  sort();

  /**
   * Write or read the data of this object to or from a stream for the purpose
   * of serialization.
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
constexpr TableIndices<N>::TableIndices(const std::size_t index0)
  : indices{index0}
{
  static_assert(
    N == 1, "This constructor is only available for TableIndices<1> objects.");
}



template <int N>
constexpr TableIndices<N>::TableIndices(const std::size_t index0,
                                        const std::size_t index1)
  : indices{index0, index1}
{
  static_assert(
    N == 2, "This constructor is only available for TableIndices<2> objects.");
}



template <int N>
constexpr TableIndices<N>::TableIndices(const std::size_t index0,
                                        const std::size_t index1,
                                        const std::size_t index2)
  : indices{index0, index1, index2}
{
  static_assert(
    N == 3, "This constructor is only available for TableIndices<3> objects.");
}



template <int N>
constexpr TableIndices<N>::TableIndices(const std::size_t index0,
                                        const std::size_t index1,
                                        const std::size_t index2,
                                        const std::size_t index3)
  : indices{index0, index1, index2, index3}
{
  static_assert(
    N == 4, "This constructor is only available for TableIndices<4> objects.");
}



template <int N>
constexpr TableIndices<N>::TableIndices(const std::size_t index0,
                                        const std::size_t index1,
                                        const std::size_t index2,
                                        const std::size_t index3,
                                        const std::size_t index4)
  : indices{index0, index1, index2, index3, index4}
{
  static_assert(
    N == 5, "This constructor is only available for TableIndices<5> objects.");
}



template <int N>
TableIndices<N>::TableIndices(const std::size_t index0,
                              const std::size_t index1,
                              const std::size_t index2,
                              const std::size_t index3,
                              const std::size_t index4,
                              const std::size_t index5,
                              const std::size_t index6,
                              const std::size_t index7,
                              const std::size_t index8)
{
  switch (N)
    {
      case 1:
        Assert(index1 == numbers::invalid_unsigned_int,
               ExcMessage("more than N index values provided"));
        DEAL_II_FALLTHROUGH;
      case 2:
        Assert(index2 == numbers::invalid_unsigned_int,
               ExcMessage("more than N index values provided"));
        DEAL_II_FALLTHROUGH;
      case 3:
        Assert(index3 == numbers::invalid_unsigned_int,
               ExcMessage("more than N index values provided"));
        DEAL_II_FALLTHROUGH;
      case 4:
        Assert(index4 == numbers::invalid_unsigned_int,
               ExcMessage("more than N index values provided"));
        DEAL_II_FALLTHROUGH;
      case 5:
        Assert(index5 == numbers::invalid_unsigned_int,
               ExcMessage("more than N index values provided"));
        DEAL_II_FALLTHROUGH;
      case 6:
        Assert(index6 == numbers::invalid_unsigned_int,
               ExcMessage("more than N index values provided"));
        DEAL_II_FALLTHROUGH;
      case 7:
        Assert(index7 == numbers::invalid_unsigned_int,
               ExcMessage("more than N index values provided"));
        DEAL_II_FALLTHROUGH;
      case 8:
        Assert(index8 == numbers::invalid_unsigned_int,
               ExcMessage("more than N index values provided"));
        break;
      default:;
    }

  // Always access "indices" with indices modulo N to avoid bogus compiler
  // warnings (although such access is always in dead code...
  switch (N)
    {
      default:
        // For TableIndices of size 10 or larger als default initialize the
        // remaining indices to numbers::invalid_unsigned_int:
        for (unsigned int i = 0; i < N; ++i)
          indices[i] = numbers::invalid_unsigned_int;
        DEAL_II_FALLTHROUGH;
      case 9:
        indices[8 % N] = index8;
        DEAL_II_FALLTHROUGH;
      case 8:
        indices[7 % N] = index7;
        DEAL_II_FALLTHROUGH;
      case 7:
        indices[6 % N] = index6;
        DEAL_II_FALLTHROUGH;
      case 6:
        indices[5 % N] = index5;
        DEAL_II_FALLTHROUGH;
      case 5:
        indices[4 % N] = index4;
        DEAL_II_FALLTHROUGH;
      case 4:
        indices[3 % N] = index3;
        DEAL_II_FALLTHROUGH;
      case 3:
        indices[2 % N] = index2;
        DEAL_II_FALLTHROUGH;
      case 2:
        indices[1 % N] = index1;
        DEAL_II_FALLTHROUGH;
      case 1:
        indices[0 % N] = index0;
    }
}


template <int N>
DEAL_II_CONSTEXPR inline std::size_t TableIndices<N>::
                                     operator[](const unsigned int i) const
{
  AssertIndexRange(i, N);
  return indices[i];
}


template <int N>
DEAL_II_CONSTEXPR inline std::size_t &TableIndices<N>::
                                      operator[](const unsigned int i)
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
DEAL_II_CONSTEXPR inline void
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
