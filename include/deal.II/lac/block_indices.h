// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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

#ifndef __deal2__block_indices_h
#define __deal2__block_indices_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <cstddef>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * @brief Auxiliary class aiding in the handling of block structures like in
 * BlockVector or FESystem.
 *
 * The information obtained from this class falls into two
 * groups. First, it is possible to obtain the number of blocks,
 * namely size(), the block_size() for each block and the total_size()
 * of the object described by the block indices, namely the length of
 * the whole index set. These functions do not make any assumption on
 * the ordering of the index set.
 *
 * If on the other hand the index set is ordered "by blocks", such
 * that each block forms a consecutive set of indices, this
 * class that manages the conversion of global indices into a block vector or
 * matrix to the local indices within this block. This is required, for
 * example, when you address a global element in a block vector and want to
 * know which element within which block this is. It is also useful if a
 * matrix is composed of several blocks, where you have to translate global
 * row and column indices to local ones.
 *
 * @ingroup data
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Wolfgang Bangerth, Guido Kanschat, 2000, 2007, 2011
 */
class BlockIndices : public Subscriptor
{
public:
  /**
   * Declare the type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Default
   * constructor. Initialize for
   * zero blocks.
   */
  BlockIndices ();

  /**
   * Constructor. Initialize the
   * number of entries in each
   * block @p i as <tt>n[i]</tt>. The
   * number of blocks will be the
   * size of the vector
   */
  BlockIndices (const std::vector<size_type> &n);

  /**
   * Specialized constructor for a
   * structure with blocks of equal size.
   */
  explicit BlockIndices(const unsigned int n_blocks, const size_type block_size = 0);

  /**
   * Reinitialize the number of
   * blocks and assign each block
   * the same number of elements.
   */
  void reinit (const unsigned int n_blocks,
               const size_type n_elements_per_block);

  /**
   * Reinitialize the number of
   * indices within each block from
   * the given argument. The number
   * of blocks will be adjusted to
   * the size of @p n and the size
   * of block @p i is set to
   * <tt>n[i]</tt>.
   */
  void reinit (const std::vector<size_type> &n);

  /**
   * Add another block of given
   * size to the end of the block
   * structure.
   */
  void push_back(const size_type size);

  /**
   * @name Size information
   */
  //@{

  /**
   * Number of blocks in index field.
   */
  unsigned int size () const;

  /**
   * Return the total number of
   * indices accumulated over all
   * blocks, that is, the dimension
   * of the vector space of the
   * block vector.
   */
  size_type total_size () const;

  /**
   * The size of the @p ith block.
   */
  size_type block_size (const unsigned int i) const;

  //@}

  /**
   * @name Index conversion
   *
   * Functions in this group
   * assume an object, which
   * was created after sorting by
   * block, such that each block
   * forms a set of consecutive
   * indices in the object.
   * If applied to other objects,
   * the numbers obtained from
   * these functions are meaningless.
   */
  //@{

  /**
   * Return the block and the
   * index within that block
   * for the global index @p i. The
   * first element of the pair is
   * the block, the second the
   * index within it.
   */
  std::pair<unsigned int,size_type>
  global_to_local (const size_type i) const;

  /**
   * Return the global index of
   * @p index in block @p block.
   */
  size_type local_to_global (const unsigned int block,
                             const size_type index) const;

  /**
   * The start index of the ith block.
   */
  size_type block_start (const unsigned int i) const;
  //@}

  /**
   * Copy operator.
   */
  BlockIndices &operator = (const BlockIndices &b);

  /**
   * Compare whether two objects
   * are the same, i.e. whether the
   * number of blocks and the sizes
   * of all blocks are equal.
   */
  bool operator == (const BlockIndices &b) const;

  /**
   * Swap the contents of these two
   * objects.
   */
  void swap (BlockIndices &b);

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

private:
  /**
   * Number of blocks. While this
   * value could be obtained
   * through
   * <tt>start_indices.size()-1</tt>,
   * we cache this value for faster
   * access.
   */
  unsigned int n_blocks;

  /**
   * Global starting index of each
   * vector. The last and redundant
   * value is the total number of
   * entries.
   */
  std::vector<size_type> start_indices;
};


/**
 * Operator for logging BlockIndices. Writes the number of blocks, the
 * size of each block and the total size of the index field.
 *
 * @ref BlockIndices
 * @author Guido Kanschat
 * @date 2011
 */
inline
LogStream &
operator << (LogStream &s, const BlockIndices &bi)
{
  const unsigned int n = bi.size();
  s << n << ":[";
  // Write first size without leading space
  if (n>0)
    s << bi.block_size(0);
  // Write all other sizes
  for (unsigned int i=1; i<n; ++i)
    s << ' ' << bi.block_size(i);
  s << "]->" << bi.total_size();
  return s;
}


template <typename MatrixType> class BlockMatrixBase;
template <typename SparsityType> class BlockSparsityPatternBase;
template <typename number>     class BlockSparseMatrixEZ;

/**
 * A class that can be used to determine whether a given type is a block
 * matrix type or not. For example,
 * @code
 *   IsBlockMatrix<SparseMatrix<double> >::value
 * @endcode
 * has the value false, whereas
 * @code
 *   IsBlockMatrix<BlockSparseMatrix<double> >::value
 * @endcode
 * is true. This is sometimes useful in template contexts where we may
 * want to do things differently depending on whether a template type
 * denotes a regular or a block matrix type.
 *
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Wolfgang Bangerth, 2009
 */
template <typename MatrixType>
struct IsBlockMatrix
{
private:
  struct yes_type
  {
    char c[1];
  };
  struct no_type
  {
    char c[2];
  };

  /**
   * Overload returning true if the class
   * is derived from BlockMatrixBase,
   * which is what block matrices do
   * (with the exception of
   * BlockSparseMatrixEZ).
   */
  template <typename T>
  static yes_type check_for_block_matrix (const BlockMatrixBase<T> *);

  /**
   * Overload returning true if the class
   * is derived from
   * BlockSparsityPatternBase, which is
   * what block sparsity patterns do.
   */
  template <typename T>
  static yes_type check_for_block_matrix (const BlockSparsityPatternBase<T> *);

  /**
   * Overload for BlockSparseMatrixEZ,
   * which is the only block matrix not
   * derived from BlockMatrixBase at the
   * time of writing this class.
   */
  template <typename T>
  static yes_type check_for_block_matrix (const BlockSparseMatrixEZ<T> *);

  /**
   * Catch all for all other potential
   * matrix types that are not block
   * matrices.
   */
  static no_type check_for_block_matrix (...);

public:
  /**
   * A statically computable value that
   * indicates whether the template
   * argument to this class is a block
   * matrix (in fact whether the type is
   * derived from BlockMatrixBase<T>).
   */
  static const bool value = (sizeof(check_for_block_matrix
                                    ((MatrixType *)0))
                             ==
                             sizeof(yes_type));
};


// instantiation of the static member
template <typename MatrixType>
const bool IsBlockMatrix<MatrixType>::value;


/* ---------------------- template and inline functions ------------------- */

inline
void
BlockIndices::reinit (const unsigned int nb,
                      const size_type block_size)
{
  n_blocks = nb;
  start_indices.resize(n_blocks+1);
  for (size_type i=0; i<=n_blocks; ++i)
    start_indices[i] = i * block_size;
}



inline
void
BlockIndices::reinit (const std::vector<size_type> &n)
{
  if (start_indices.size() != n.size()+1)
    {
      n_blocks = static_cast<unsigned int>(n.size());
      start_indices.resize(n_blocks+1);
    }
  start_indices[0] = 0;
  for (size_type i=1; i<=n_blocks; ++i)
    start_indices[i] = start_indices[i-1] + n[i-1];
}


inline
BlockIndices::BlockIndices ()
  :
  n_blocks(0),
  start_indices(1, 0)
{}



inline
BlockIndices::BlockIndices (
  const unsigned int n_blocks,
  const size_type block_size)
  :
  n_blocks(n_blocks),
  start_indices(n_blocks+1)
{
  for (size_type i=0; i<=n_blocks; ++i)
    start_indices[i] = i * block_size;
}



inline
BlockIndices::BlockIndices (const std::vector<size_type> &n)
  :
  n_blocks(static_cast<unsigned int>(n.size())),
  start_indices(n.size()+1)
{
  reinit (n);
}


inline
void
BlockIndices::push_back(const size_type sz)
{
  start_indices.push_back(start_indices[n_blocks]+sz);
  ++n_blocks;
  AssertDimension(start_indices.size(), n_blocks+1);
}


inline
std::pair<unsigned int,BlockIndices::size_type>
BlockIndices::global_to_local (const size_type i) const
{
  Assert (i<total_size(), ExcIndexRangeType<size_type>(i, 0, total_size()));
  Assert (n_blocks > 0, ExcLowerRangeType<size_type>(i, size_type(1)));

  unsigned int block = n_blocks-1;
  while (i < start_indices[block])
    --block;

  return std::pair<size_type,size_type>(block,
                                        i-start_indices[block]);
}


inline
BlockIndices::size_type
BlockIndices::local_to_global (const unsigned int block,
                               const size_type index) const
{
  Assert (block < n_blocks, ExcIndexRange(block, 0, n_blocks));
  Assert (index < start_indices[block+1]-start_indices[block],
          ExcIndexRangeType<size_type> (index, 0, start_indices[block+1]-start_indices[block]));

  return start_indices[block]+index;
}


inline
unsigned int
BlockIndices::size () const
{
  return n_blocks;
}



inline
BlockIndices::size_type
BlockIndices::total_size () const
{
  if (n_blocks == 0) return 0;
  return start_indices[n_blocks];
}



inline
BlockIndices::size_type
BlockIndices::block_size (const unsigned int block) const
{
  Assert (block < n_blocks, ExcIndexRange(block, 0, n_blocks));
  return start_indices[block+1]-start_indices[block];
}



inline
BlockIndices::size_type
BlockIndices::block_start (const unsigned int block) const
{
  Assert (block < n_blocks, ExcIndexRange(block, 0, n_blocks));
  return start_indices[block];
}



inline
BlockIndices &
BlockIndices::operator = (const BlockIndices &b)
{
  start_indices = b.start_indices;
  n_blocks = b.n_blocks;
  return *this;
}



inline
bool
BlockIndices::operator == (const BlockIndices &b) const
{
  if (n_blocks != b.n_blocks)
    return false;

  for (size_type i=0; i<=n_blocks; ++i)
    if (start_indices[i] != b.start_indices[i])
      return false;

  return true;
}



inline
void
BlockIndices::swap (BlockIndices &b)
{
  Assert (n_blocks == b.n_blocks,
          ExcDimensionMismatch(n_blocks, b.n_blocks));

  for (size_type i=0; i<=n_blocks; ++i)
    std::swap (start_indices[i], b.start_indices[i]);
}



inline
std::size_t
BlockIndices::memory_consumption () const
{
  return (sizeof(*this) +
          start_indices.size() * sizeof(start_indices[0]));
}



/* ----------------- global functions ---------------------------- */


/**
 * Global function @p swap which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two objects.
 *
 * @relates BlockIndices
 * @author Wolfgang Bangerth, 2000
 */
inline
void swap (BlockIndices &u, BlockIndices &v)
{
  u.swap (v);
}




DEAL_II_NAMESPACE_CLOSE

#endif
