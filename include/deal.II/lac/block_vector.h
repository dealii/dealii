// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__block_vector_h
#define __deal2__block_vector_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/block_vector_base.h>

#include <cstdio>
#include <vector>

DEAL_II_NAMESPACE_OPEN


#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  class Vector;
  class BlockVector;
}
#endif



/*! @addtogroup Vectors
 *@{
 */


/**
 * An implementation of block vectors based on deal.II vectors. While the base
 * class provides for most of the interface, this class handles the actual
 * allocation of vectors and provides functions that are specific to the
 * underlying vector type.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
 *
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2000, 2001, 2002, 2004
 */
template <typename Number>
class BlockVector : public BlockVectorBase<Vector<Number> >
{
public:
  /**
   * Typedef the base class for simpler
   * access to its own typedefs.
   */
  typedef BlockVectorBase<Vector<Number> > BaseClass;

  /**
   * Typedef the type of the underlying
   * vector.
   */
  typedef typename BaseClass::BlockType  BlockType;

  /**
   * Import the typedefs from the base
   * class.
   */
  typedef typename BaseClass::value_type      value_type;
  typedef typename BaseClass::real_type       real_type;
  typedef typename BaseClass::pointer         pointer;
  typedef typename BaseClass::const_pointer   const_pointer;
  typedef typename BaseClass::reference       reference;
  typedef typename BaseClass::const_reference const_reference;
  typedef typename BaseClass::size_type       size_type;
  typedef typename BaseClass::iterator        iterator;
  typedef typename BaseClass::const_iterator  const_iterator;

  /**
   *  Constructor. There are three
   *  ways to use this
   *  constructor. First, without
   *  any arguments, it generates
   *  an object with no
   *  blocks. Given one argument,
   *  it initializes <tt>num_blocks</tt>
   *  blocks, but these blocks have
   *  size zero. The third variant
   *  finally initializes all
   *  blocks to the same size
   *  <tt>block_size</tt>.
   *
   *  Confer the other constructor
   *  further down if you intend to
   *  use blocks of different
   *  sizes.
   */
  explicit BlockVector (const unsigned int num_blocks = 0,
                        const size_type block_size = 0);

  /**
   * Copy-Constructor. Dimension set to
   * that of V, all components are copied
   * from V
   */
  BlockVector (const BlockVector<Number> &V);


#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG
  /**
   * Copy constructor taking a BlockVector of
   * another data type. This will fail if
   * there is no conversion path from
   * <tt>OtherNumber</tt> to <tt>Number</tt>. Note that
   * you may lose accuracy when copying
   * to a BlockVector with data elements with
   * less accuracy.
   *
   * Older versions of gcc did not honor
   * the @p explicit keyword on template
   * constructors. In such cases, it is
   * easy to accidentally write code that
   * can be very inefficient, since the
   * compiler starts performing hidden
   * conversions. To avoid this, this
   * function is disabled if we have
   * detected a broken compiler during
   * configuration.
   */
  template <typename OtherNumber>
  explicit
  BlockVector (const BlockVector<OtherNumber> &v);
#endif


#ifdef DEAL_II_WITH_TRILINOS
  /**
   * A copy constructor taking a
   * (parallel) Trilinos block
   * vector and copying it into
   * the deal.II own format.
   */
  BlockVector (const TrilinosWrappers::BlockVector &v);

#endif
  /**
   * Constructor. Set the number of
   * blocks to
   * <tt>block_sizes.size()</tt> and
   * initialize each block with
   * <tt>block_sizes[i]</tt> zero
   * elements.
   */
  BlockVector (const std::vector<size_type> &block_sizes);

  /**
   * Constructor. Initialize vector
   * to the structure found in the
   * BlockIndices argument.
   */
  BlockVector (const BlockIndices &block_indices);

  /**
   * Constructor. Set the number of
   * blocks to
   * <tt>n.size()</tt>. Initialize the
   * vector with the elements
   * pointed to by the range of
   * iterators given as second and
   * third argument. Apart from the
   * first argument, this
   * constructor is in complete
   * analogy to the respective
   * constructor of the
   * <tt>std::vector</tt> class, but the
   * first argument is needed in
   * order to know how to subdivide
   * the block vector into
   * different blocks.
   */
  template <typename InputIterator>
  BlockVector (const std::vector<size_type>    &n,
               const InputIterator              first,
               const InputIterator              end);

  /**
   * Destructor. Clears memory
   */
  ~BlockVector ();

  /**
   * Call the compress() function on all
   * the subblocks.
  *
  * This functionality only needs to be
  * called if using MPI based vectors and
  * exists in other objects for
  * compatibility.
  *
  * See @ref GlossCompress "Compressing
  * distributed objects" for more
  * information.
   */
  void compress (::dealii::VectorOperation::values operation
                 =::dealii::VectorOperation::unknown);

  /**
   * Copy operator: fill all components of
   * the vector with the given scalar
   * value.
   */
  BlockVector &operator = (const value_type s);

  /**
   * Copy operator for arguments of the
   * same type. Resize the
   * present vector if necessary.
   */
  BlockVector &
  operator= (const BlockVector &V);

  /**
   * Copy operator for template arguments
   * of different types. Resize the
   * present vector if necessary.
   */
  template <class Number2>
  BlockVector &
  operator= (const BlockVector<Number2> &V);

  /**
   * Copy a regular vector into a
   * block vector.
   */
  BlockVector &
  operator= (const Vector<Number> &V);

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * A copy constructor from a
   * Trilinos block vector to a
   * deal.II block vector.
   */
  BlockVector &
  operator= (const TrilinosWrappers::BlockVector &V);
#endif
  /**
   * Reinitialize the BlockVector to
   * contain <tt>num_blocks</tt> blocks of
   * size <tt>block_size</tt> each.
   *
   * If the second argument is left
   * at its default value, then the
   * block vector allocates the
   * specified number of blocks but
   * leaves them at zero size. You
   * then need to later
   * reinitialize the individual
   * blocks, and call
   * collect_sizes() to update the
   * block system's knowledge of
   * its individual block's sizes.
   *
   * If <tt>fast==false</tt>, the vector
   * is filled with zeros.
   */
  void reinit (const unsigned int num_blocks,
               const size_type block_size = 0,
               const bool fast = false);

  /**
   * Reinitialize the BlockVector such that
   * it contains
   * <tt>block_sizes.size()</tt>
   * blocks. Each block is reinitialized to
   * dimension <tt>block_sizes[i]</tt>.
   *
   * If the number of blocks is the
   * same as before this function
   * was called, all vectors remain
   * the same and reinit() is
   * called for each vector.
   *
   * If <tt>fast==false</tt>, the vector
   * is filled with zeros.
   *
   * Note that you must call this
   * (or the other reinit()
   * functions) function, rather
   * than calling the reinit()
   * functions of an individual
   * block, to allow the block
   * vector to update its caches of
   * vector sizes. If you call
   * reinit() on one of the
   * blocks, then subsequent
   * actions on this object may
   * yield unpredictable results
   * since they may be routed to
   * the wrong block.
   */
  void reinit (const std::vector<size_type> &N,
               const bool                    fast=false);

  /**
   * Reinitialize the BlockVector
   * to reflect the structure found
   * in BlockIndices.
   *
   * If the number of blocks is the
   * same as before this function
   * was called, all vectors remain
   * the same and reinit() is
   * called for each vector.
   *
   * If <tt>fast==false</tt>, the vector
   * is filled with zeros.
   */
  void reinit (const BlockIndices &block_indices,
               const bool fast=false);

  /**
   * Change the dimension to that
   * of the vector <tt>V</tt>. The same
   * applies as for the other
   * reinit() function.
   *
   * The elements of <tt>V</tt> are not
   * copied, i.e.  this function is
   * the same as calling <tt>reinit
   * (V.size(), fast)</tt>.
   *
   * Note that you must call this
   * (or the other reinit()
   * functions) function, rather
   * than calling the reinit()
   * functions of an individual
   * block, to allow the block
   * vector to update its caches of
   * vector sizes. If you call
   * reinit() of one of the
   * blocks, then subsequent
   * actions of this object may
   * yield unpredictable results
   * since they may be routed to
   * the wrong block.
   */
  template <typename Number2>
  void reinit (const BlockVector<Number2> &V,
               const bool                 fast=false);

  /**
   * Scale each element of the
   * vector by the given factor.
   *
   * This function is deprecated
   * and will be removed in a
   * future version. Use
   * <tt>operator *=</tt> and
   * <tt>operator /=</tt> instead.
   *
   * @deprecated Use <tt>operator*=</tt>
   * instead.
   */
  void scale (const value_type factor) DEAL_II_DEPRECATED;

  /**
   * Multiply each element of this
   * vector by the corresponding
   * element of <tt>v</tt>.
   */
  template <class BlockVector2>
  void scale (const BlockVector2 &v);

  /**
   * Swap the contents of this
   * vector and the other vector
   * <tt>v</tt>. One could do this
   * operation with a temporary
   * variable and copying over the
   * data elements, but this
   * function is significantly more
   * efficient since it only swaps
   * the pointers to the data of
   * the two vectors and therefore
   * does not need to allocate
   * temporary storage and move
   * data around.
   *
   * Limitation: right now this
   * function only works if both
   * vectors have the same number
   * of blocks. If needed, the
   * numbers of blocks should be
   * exchanged, too.
   *
   * This function is analog to the
   * the swap() function of all C++
   * standard containers. Also,
   * there is a global function
   * swap(u,v) that simply calls
   * <tt>u.swap(v)</tt>, again in analogy
   * to standard functions.
   */
  void swap (BlockVector<Number> &v);

  /**
   *  Output of vector in user-defined
   *  format.
   */
  void print (const char *format = 0) const;

  /**
   * Print to a stream.
   */
  void print (std::ostream       &out,
              const unsigned int  precision = 3,
              const bool          scientific = true,
              const bool          across = true) const;

  /**
   * Write the vector en bloc to a
   * stream. This is done in a binary mode,
   * so the output is neither readable by
   * humans nor (probably) by other
   * computers using a different operating
   * system or number format.
   */
  void block_write (std::ostream &out) const;

  /**
   * Read a vector en block from a
   * file. This is done using the inverse
   * operations to the above function, so
   * it is reasonably fast because the
   * bitstream is not interpreted.
   *
   * The vector is resized if necessary.
   *
   * A primitive form of error checking is
   * performed which will recognize the
   * bluntest attempts to interpret some
   * data as a vector stored bitwise to a
   * file, but not more.
   */
  void block_read (std::istream &in);

  /** @addtogroup Exceptions
   * @{ */

  /**
   * Exception
   */
  DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);
  //@}
};

/*@}*/

#ifndef DOXYGEN
/*----------------------- Inline functions ----------------------------------*/



template <typename Number>
template <typename InputIterator>
BlockVector<Number>::BlockVector (const std::vector<size_type>    &n,
                                  const InputIterator              first,
                                  const InputIterator              end)
{
  // first set sizes of blocks, but
  // don't initialize them as we will
  // copy elements soon
  reinit (n, true);
  InputIterator start = first;
  for (size_type b=0; b<n.size(); ++b)
    {
      InputIterator end = start;
      std::advance (end, static_cast<signed int>(n[b]));
      std::copy (start, end, this->block(b).begin());
      start = end;
    };
  Assert (start == end, ExcIteratorRangeDoesNotMatchVectorSize());
}



template <typename Number>
inline
BlockVector<Number> &
BlockVector<Number>::operator = (const value_type s)
{

  Assert (numbers::is_finite(s), ExcNumberNotFinite());

  BaseClass::operator = (s);
  return *this;
}



template <typename Number>
inline
BlockVector<Number> &
BlockVector<Number>::operator = (const BlockVector &v)
{
  reinit (v, true);
  BaseClass::operator = (v);
  return *this;
}



template <typename Number>
inline
BlockVector<Number> &
BlockVector<Number>::operator = (const Vector<Number> &v)
{
  BaseClass::operator = (v);
  return *this;
}



template <typename Number>
template <typename Number2>
inline
BlockVector<Number> &
BlockVector<Number>::operator = (const BlockVector<Number2> &v)
{
  reinit (v, true);
  BaseClass::operator = (v);
  return *this;
}

template <typename Number>
inline
void BlockVector<Number>::compress (::dealii::VectorOperation::values operation)
{
  for (size_type i=0; i<this->n_blocks(); ++i)
    this->components[i].compress(operation);
}



template <typename Number>
void BlockVector<Number>::scale (const value_type factor)
{

  Assert (numbers::is_finite(factor), ExcNumberNotFinite());

  for (size_type i=0; i<this->n_blocks(); ++i)
    this->components[i] *= factor;
}



template <typename Number>
template <class BlockVector2>
void BlockVector<Number>::scale (const BlockVector2 &v)
{
  BaseClass::scale (v);
}

#endif // DOXYGEN


/**
 * Global function which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @relates BlockVector
 * @author Wolfgang Bangerth, 2000
 */
template <typename Number>
inline
void swap (BlockVector<Number> &u,
           BlockVector<Number> &v)
{
  u.swap (v);
}

DEAL_II_NAMESPACE_CLOSE

#endif
