// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2016 by the deal.II authors
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

#ifndef dealii__table_indices_h
#define dealii__table_indices_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <algorithm>
#include <ostream>
#include <iterator>


DEAL_II_NAMESPACE_OPEN


/**
 * A class representing a fixed size array of indices.
 *
 * It is used in tensorial objects like the TableBase and SymmetricTensor
 * classes to represent a nested choice of indices.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, Matthias Maier, 2002, 2015
 */
template <int N>
class TableIndices
{
public:

  /**
   * Default constructor. It sets all indices to zero.
   */
  TableIndices();

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
   * @tparam N The number of indices stored in each object.
   */
  TableIndices (const unsigned int index0,
                const unsigned int index1 = numbers::invalid_unsigned_int,
                const unsigned int index2 = numbers::invalid_unsigned_int,
                const unsigned int index3 = numbers::invalid_unsigned_int,
                const unsigned int index4 = numbers::invalid_unsigned_int,
                const unsigned int index5 = numbers::invalid_unsigned_int,
                const unsigned int index6 = numbers::invalid_unsigned_int,
                const unsigned int index7 = numbers::invalid_unsigned_int,
                const unsigned int index8 = numbers::invalid_unsigned_int);

  /**
   * Read-only access the value of the <tt>i</tt>th index.
   */
  unsigned int operator[] (const unsigned int i) const;

  /**
   * Write access the value of the <tt>i</tt>th index.
   */
  unsigned int &operator[] (const unsigned int i);

  /**
   * Compare two index fields for equality.
   */
  bool operator == (const TableIndices<N> &other) const;

  /**
   * Compare two index fields for inequality.
   */
  bool operator != (const TableIndices<N> &other) const;

  /**
   * Sort the indices in ascending order. While this operation is not very
   * useful for Table objects, it is used for the SymmetricTensor class.
   */
  void sort ();

  /**
   * Write or read the data of this object to or from a stream for the purpose
   * of serialization.
   */
  template <class Archive>
  void serialize (Archive &ar, const unsigned int version);

protected:
  /**
   * Store the indices in an array.
   */
  unsigned int indices[N];
};



/* --------------------- Template and inline functions ---------------- */


template <int N>
TableIndices<N>::TableIndices()
{
  Assert (N > 0, ExcMessage("Cannot create a TableIndices object of size 0"));

  for (unsigned int i=0; i<N; ++i)
    indices[i] = 0;
}


template <int N>
TableIndices<N>::TableIndices(const unsigned int index0,
                              const unsigned int index1,
                              const unsigned int index2,
                              const unsigned int index3,
                              const unsigned int index4,
                              const unsigned int index5,
                              const unsigned int index6,
                              const unsigned int index7,
                              const unsigned int index8)
{
  Assert (N > 0, ExcMessage("Cannot create a TableIndices object of size 0"));

  switch (N)
    {
    case 1:
      Assert (index1 == numbers::invalid_unsigned_int, ExcMessage("more than N index values provided"));
      DEAL_II_FALLTHROUGH;
    case 2:
      Assert (index2 == numbers::invalid_unsigned_int, ExcMessage("more than N index values provided"));
      DEAL_II_FALLTHROUGH;
    case 3:
      Assert (index3 == numbers::invalid_unsigned_int, ExcMessage("more than N index values provided"));
      DEAL_II_FALLTHROUGH;
    case 4:
      Assert (index4 == numbers::invalid_unsigned_int, ExcMessage("more than N index values provided"));
      DEAL_II_FALLTHROUGH;
    case 5:
      Assert (index5 == numbers::invalid_unsigned_int, ExcMessage("more than N index values provided"));
      DEAL_II_FALLTHROUGH;
    case 6:
      Assert (index6 == numbers::invalid_unsigned_int, ExcMessage("more than N index values provided"));
      DEAL_II_FALLTHROUGH;
    case 7:
      Assert (index7 == numbers::invalid_unsigned_int, ExcMessage("more than N index values provided"));
      DEAL_II_FALLTHROUGH;
    case 8:
      Assert (index8 == numbers::invalid_unsigned_int, ExcMessage("more than N index values provided"));
      break;
    default:
      ;
    }

  // Always access "indices" with indices modulo N to avoid bogus compiler
  // warnings (although such access is always in dead code...
  switch (N)
    {
    default:
      // For TableIndices of size 10 or larger als default initialize the
      // remaining indices to numbers::invalid_unsigned_int:
      for (unsigned int i=0; i<N; ++i)
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
inline
unsigned int
TableIndices<N>::operator [] (const unsigned int i) const
{
  Assert (i < N, ExcIndexRange (i, 0, N));
  return indices[i];
}


template <int N>
inline
unsigned int &
TableIndices<N>::operator [] (const unsigned int i)
{
  Assert (i < N, ExcIndexRange (i, 0, N));
  return indices[i];
}


template <int N>
inline
bool
TableIndices<N>::operator == (const TableIndices<N> &other) const
{
  for (unsigned int i=0; i<N; ++i)
    if (indices[i] != other.indices[i])
      return false;
  return true;
}


template <int N>
inline
bool
TableIndices<N>::operator != (const TableIndices<N> &other) const
{
  return !(*this == other);
}


template <int N>
inline
void
TableIndices<N>::sort ()
{
  std::sort(std::begin(indices), std::end(indices));
}


template <int N>
template <class Archive>
inline
void
TableIndices<N>::serialize (Archive &ar, const unsigned int)
{
  ar &indices;
}


/**
 * Output operator for TableIndices objects; reports them in a list like this:
 * <code>[i1,i2,...]</code>.
 *
 * @relates TableIndices
 */
template <int N>
std::ostream &
operator << (std::ostream &out,
             const TableIndices<N> &indices)
{
  out << '[';
  for (unsigned int i=0; i<N; ++i)
    {
      out << indices[i];
      if (i+1 != N)
        out << ',';
    }
  out << ']';

  return out;
}


DEAL_II_NAMESPACE_CLOSE

#endif
