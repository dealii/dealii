// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2002 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_block_sparse_matrix_ez_h
#define dealii_block_sparse_matrix_ez_h


// TODO: Derive BlockSparseMatrixEZ from BlockMatrixBase, like all the
// other block matrices as well; this would allow to instantiate a few
// functions with this template argument as well (in particular
// AffineConstraints::distribute_local_to_global)

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/observer_pointer.h>
#include <deal.II/base/table.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/block_matrix_base.h>
#include <deal.II/lac/sparse_matrix_ez.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <typename Number>
class BlockVector;
#endif

/**
 * @addtogroup Matrix1
 * @{
 */


/**
 * A block matrix consisting of blocks of type SparseMatrixEZ.
 *
 * Like the other Block-objects, this matrix can be used like a
 * SparseMatrixEZ, when it comes to access to entries. Then, there are
 * functions for the multiplication with BlockVector and access to the
 * individual blocks.
 *
 * @see
 * @ref GlossBlockLA "Block (linear algebra)"
 */
template <typename Number>
class BlockSparseMatrixEZ : public BlockMatrixBase<SparseMatrixEZ<Number>>
{
public:
  /**
   * Typedef the base class for simpler access to its own alias.
   */
  using BaseClass = BlockMatrixBase<SparseMatrixEZ<Number>>;

  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Default constructor. The result is an empty object with zero dimensions.
   */
  BlockSparseMatrixEZ() = default;

  /**
   * Constructor setting up an object with given number of block rows and
   * columns. The blocks themselves still have zero dimension.
   */
  BlockSparseMatrixEZ(const unsigned int block_rows,
                      const unsigned int block_cols);

  /**
   * Copy constructor. This is needed for some container classes. It creates
   * an object of the same number of block rows and columns. To avoid unwanted
   * loss of information, the blocks of the other matrix must be empty.
   */
  BlockSparseMatrixEZ(const BlockSparseMatrixEZ<Number> &);

  /**
   * Copy operator. Like the copy constructor, this may be called for objects
   * with empty blocks only.
   */
  BlockSparseMatrixEZ &
  operator=(const BlockSparseMatrixEZ<Number> &);

  /**
   * This operator assigns a scalar to a matrix. Since this does usually not
   * make much sense (should we set all matrix entries to this value? Only the
   * nonzero entries of the sparsity pattern?), this operation is only allowed
   * if the actual value to be assigned is zero. This operator only exists to
   * allow for the obvious notation <tt>matrix=0</tt>, which sets all elements
   * of the matrix to zero, but keep the sparsity pattern previously used.
   */
  BlockSparseMatrixEZ &
  operator=(const double d);

  /**
   * Set matrix to zero dimensions and release memory.
   */
  void
  clear();

  /**
   * Initialize to given block numbers.  After this operation, the matrix will
   * have the block dimensions provided. Each block will have zero dimensions
   * and must be initialized subsequently. After setting the sizes of the
   * blocks, collect_sizes() must be called to update internal data
   * structures.
   */
  void
  reinit(const unsigned int n_block_rows, const unsigned int n_block_cols);

  /**
   * This function collects the sizes of the sub-objects and stores them in
   * internal arrays, in order to be able to relay global indices into the
   * matrix to indices into the subobjects. You *must* call this function each
   * time after you have changed the size of the sub-objects.
   */
  void
  collect_sizes();

  /**
   * Return whether the object is empty. It is empty if no memory is
   * allocated, which is the same as that both dimensions are zero. This
   * function is just the concatenation of the respective call to all
   * sub-matrices.
   */
  bool
  empty() const;

  /**
   * @name Multiplications
   */
  /** @{ */
  /**
   * Matrix-vector multiplication: let $dst = M*src$ with $M$ being this
   * matrix.
   */
  template <typename block_number>
  void
  vmult(BlockVector<block_number>       &dst,
        const BlockVector<block_number> &src) const;

  /**
   * Matrix-vector multiplication. Just like the previous function, but only
   * applicable if the matrix has only one block column.
   */
  template <typename block_number, typename nonblock_number>
  void
  vmult(BlockVector<block_number>     &dst,
        const Vector<nonblock_number> &src) const;

  /**
   * Matrix-vector multiplication. Just like the previous function, but only
   * applicable if the matrix has only one block row.
   */
  template <typename block_number, typename nonblock_number>
  void
  vmult(Vector<nonblock_number>         &dst,
        const BlockVector<block_number> &src) const;

  /**
   * Matrix-vector multiplication. Just like the previous function, but only
   * applicable if the matrix has only one block.
   */
  template <typename nonblock_number>
  void
  vmult(Vector<nonblock_number> &dst, const Vector<nonblock_number> &src) const;

  /**
   * Matrix-vector multiplication: let $dst = M^T*src$ with $M$ being this
   * matrix. This function does the same as vmult() but takes the transposed
   * matrix.
   */
  template <typename block_number>
  void
  Tvmult(BlockVector<block_number>       &dst,
         const BlockVector<block_number> &src) const;

  /**
   * Matrix-vector multiplication. Just like the previous function, but only
   * applicable if the matrix has only one block row.
   */
  template <typename block_number, typename nonblock_number>
  void
  Tvmult(BlockVector<block_number>     &dst,
         const Vector<nonblock_number> &src) const;

  /**
   * Matrix-vector multiplication. Just like the previous function, but only
   * applicable if the matrix has only one block column.
   */
  template <typename block_number, typename nonblock_number>
  void
  Tvmult(Vector<nonblock_number>         &dst,
         const BlockVector<block_number> &src) const;

  /**
   * Matrix-vector multiplication. Just like the previous function, but only
   * applicable if the matrix has only one block.
   */
  template <typename nonblock_number>
  void
  Tvmult(Vector<nonblock_number>       &dst,
         const Vector<nonblock_number> &src) const;
  /** @} */

  /**
   * Print statistics. If @p full is @p true, prints a histogram of all
   * existing row lengths and allocated row lengths. Otherwise, just the
   * relation of allocated and used entries is shown.
   */
  template <typename StreamType>
  void
  print_statistics(StreamType &s, bool full = false);
};

/** @} */
/*----------------------------------------------------------------------*/



template <typename number>
template <typename block_number>
inline void
BlockSparseMatrixEZ<number>::vmult(BlockVector<block_number>       &dst,
                                   const BlockVector<block_number> &src) const
{
  BaseClass::vmult_block_block(dst, src);
}

template <typename number>
template <typename block_number, typename nonblock_number>
inline void
BlockSparseMatrixEZ<number>::vmult(BlockVector<block_number>     &dst,
                                   const Vector<nonblock_number> &src) const
{
  BaseClass::vmult_block_nonblock(dst, src);
}

template <typename number>
template <typename block_number, typename nonblock_number>
inline void
BlockSparseMatrixEZ<number>::vmult(Vector<nonblock_number>         &dst,
                                   const BlockVector<block_number> &src) const
{
  BaseClass::vmult_nonblock_block(dst, src);
}

template <typename number>
template <typename nonblock_number>
inline void
BlockSparseMatrixEZ<number>::vmult(Vector<nonblock_number>       &dst,
                                   const Vector<nonblock_number> &src) const
{
  BaseClass::vmult_nonblock_nonblock(dst, src);
}

template <typename number>
template <typename block_number>
inline void
BlockSparseMatrixEZ<number>::Tvmult(BlockVector<block_number>       &dst,
                                    const BlockVector<block_number> &src) const
{
  BaseClass::Tvmult_block_block(dst, src);
}

template <typename number>
template <typename block_number, typename nonblock_number>
inline void
BlockSparseMatrixEZ<number>::Tvmult(BlockVector<block_number>     &dst,
                                    const Vector<nonblock_number> &src) const
{
  BaseClass::Tvmult_block_nonblock(dst, src);
}

template <typename number>
template <typename block_number, typename nonblock_number>
inline void
BlockSparseMatrixEZ<number>::Tvmult(Vector<nonblock_number>         &dst,
                                    const BlockVector<block_number> &src) const
{
  BaseClass::Tvmult_nonblock_block(dst, src);
}

template <typename number>
template <typename nonblock_number>
inline void
BlockSparseMatrixEZ<number>::Tvmult(Vector<nonblock_number>       &dst,
                                    const Vector<nonblock_number> &src) const
{
  BaseClass::Tvmult_nonblock_nonblock(dst, src);
}



template <typename number>
template <typename StreamType>
inline void
BlockSparseMatrixEZ<number>::print_statistics(StreamType &out, bool full)
{
  size_type              used_total      = 0;
  size_type              allocated_total = 0;
  size_type              reserved_total  = 0;
  std::vector<size_type> used_by_line_total;

  size_type              used;
  size_type              allocated;
  size_type              reserved;
  std::vector<size_type> used_by_line;

  for (size_type i = 0; i < this->n_block_rows(); ++i)
    for (size_type j = 0; j < this->n_block_cols(); ++j)
      {
        used_by_line.clear();
        out << "block:\t" << i << '\t' << j << std::endl;
        this->sub_objects(i, j)->compute_statistics(
          used, allocated, reserved, used_by_line, full);

        out << "used:" << used << std::endl
            << "allocated:" << allocated << std::endl
            << "reserved:" << reserved << std::endl;

        used_total += used;
        allocated_total += allocated;
        reserved_total += reserved;

        if (full)
          {
            used_by_line_total.resize(used_by_line.size());
            for (size_type i = 0; i < used_by_line.size(); ++i)
              if (used_by_line[i] != 0)
                {
                  out << "row-entries\t" << i << "\trows\t" << used_by_line[i]
                      << std::endl;
                  used_by_line_total[i] += used_by_line[i];
                }
          }
      }
  out << "Total" << std::endl
      << "used:" << used_total << std::endl
      << "allocated:" << allocated_total << std::endl
      << "reserved:" << reserved_total << std::endl;
  for (size_type i = 0; i < used_by_line_total.size(); ++i)
    if (used_by_line_total[i] != 0)
      {
        out << "row-entries\t" << i << "\trows\t" << used_by_line_total[i]
            << std::endl;
      }
}


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_block_sparse_matrix_ez_h
