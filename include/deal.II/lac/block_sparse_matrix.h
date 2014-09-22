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

#ifndef __deal2__block_sparse_matrix_h
#define __deal2__block_sparse_matrix_h


#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/lac/block_matrix_base.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/exceptions.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


/*! @addtogroup Matrix1
 *@{
 */


/**
 * Blocked sparse matrix based on the SparseMatrix class. This class
 * implements the functions that are specific to the SparseMatrix base objects
 * for a blocked sparse matrix, and leaves the actual work relaying most of
 * the calls to the individual blocks to the functions implemented in the base
 * class. See there also for a description of when this class is useful.
 *
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Wolfgang Bangerth, 2000, 2004
 */
template <typename number>
class BlockSparseMatrix : public BlockMatrixBase<SparseMatrix<number> >
{
public:
  /**
   * Typedef the base class for simpler
   * access to its own typedefs.
   */
  typedef BlockMatrixBase<SparseMatrix<number> > BaseClass;

  /**
   * Typedef the type of the underlying
   * matrix.
   */
  typedef typename BaseClass::BlockType  BlockType;

  /**
   * Import the typedefs from the base
   * class.
   */
  typedef typename BaseClass::value_type      value_type;
  typedef typename BaseClass::pointer         pointer;
  typedef typename BaseClass::const_pointer   const_pointer;
  typedef typename BaseClass::reference       reference;
  typedef typename BaseClass::const_reference const_reference;
  typedef typename BaseClass::size_type       size_type;
  typedef typename BaseClass::iterator        iterator;
  typedef typename BaseClass::const_iterator  const_iterator;

  /**
   * @name Constructors and initalization
   */
//@{
  /**
   * Constructor; initializes the
   * matrix to be empty, without
   * any structure, i.e.  the
   * matrix is not usable at
   * all. This constructor is
   * therefore only useful for
   * matrices which are members of
   * a class. All other matrices
   * should be created at a point
   * in the data flow where all
   * necessary information is
   * available.
   *
   * You have to initialize the
   * matrix before usage with
   * reinit(BlockSparsityPattern). The
   * number of blocks per row and
   * column are then determined by
   * that function.
   */
  BlockSparseMatrix ();

  /**
   * Constructor. Takes the given
   * matrix sparsity structure to
   * represent the sparsity pattern
   * of this matrix. You can change
   * the sparsity pattern later on
   * by calling the reinit()
   * function.
   *
   * This constructor initializes
   * all sub-matrices with the
   * sub-sparsity pattern within
   * the argument.
   *
   * You have to make sure that the
   * lifetime of the sparsity
   * structure is at least as long
   * as that of this matrix or as
   * long as reinit() is not called
   * with a new sparsity structure.
   */
  BlockSparseMatrix (const BlockSparsityPattern &sparsity);

  /**
   * Destructor.
   */
  virtual ~BlockSparseMatrix ();



  /**
   * Pseudo copy operator only copying
   * empty objects. The sizes of the block
   * matrices need to be the same.
   */
  BlockSparseMatrix &
  operator = (const BlockSparseMatrix &);

  /**
   * This operator assigns a scalar to a
   * matrix. Since this does usually not
   * make much sense (should we set all
   * matrix entries to this value? Only
   * the nonzero entries of the sparsity
   * pattern?), this operation is only
   * allowed if the actual value to be
   * assigned is zero. This operator only
   * exists to allow for the obvious
   * notation <tt>matrix=0</tt>, which
   * sets all elements of the matrix to
   * zero, but keep the sparsity pattern
   * previously used.
   */
  BlockSparseMatrix &
  operator = (const double d);

  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor. It also forgets
   * the sparsity pattern it was
   * previously tied to.
   *
   * This calls SparseMatrix::clear on all
   * sub-matrices and then resets this
   * object to have no blocks at all.
   */
  void clear ();

  /**
   * Reinitialize the sparse matrix
   * with the given sparsity
   * pattern. The latter tells the
   * matrix how many nonzero
   * elements there need to be
   * reserved.
   *
   * Basically, this function only
   * calls SparseMatrix::reinit() of the
   * sub-matrices with the block
   * sparsity patterns of the
   * parameter.
   *
   * You have to make sure that the lifetime of the sparsity structure is at
   * least as long as that of this matrix or as long as reinit(const
   * SparsityPattern &) is not called with a new sparsity structure.
   *
   * The elements of the matrix are set to zero by this function.
   */
  virtual void reinit (const BlockSparsityPattern &sparsity);
//@}

  /**
   * @name Information on the matrix
   */
//@{
  /**
   * Return whether the object is
   * empty. It is empty if either
   * both dimensions are zero or no
   * BlockSparsityPattern is
   * associated.
   */
  bool empty () const;

  /**
   * Return the number of entries
   * in a specific row.
   */
  size_type get_row_length (const size_type row) const;

  /**
   * Return the number of nonzero
   * elements of this
   * matrix. Actually, it returns
   * the number of entries in the
   * sparsity pattern; if any of
   * the entries should happen to
   * be zero, it is counted anyway.
   */
  size_type n_nonzero_elements () const;

  /**
   * Return the number of actually
   * nonzero elements. Just counts the
   * number of actually nonzero elements
   * (with absolute value larger than
   * threshold) of all the blocks.
   */
  size_type n_actually_nonzero_elements (const double threshold = 0.0) const;

  /**
   * Return a (constant) reference
   * to the underlying sparsity
   * pattern of this matrix.
   *
   * Though the return value is
   * declared <tt>const</tt>, you
   * should be aware that it may
   * change if you call any
   * nonconstant function of
   * objects which operate on it.
   */
  const BlockSparsityPattern &
  get_sparsity_pattern () const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;
//@}

  /**
   * @name Multiplications
   */
//@{
  /**
   * Matrix-vector multiplication:
   * let $dst = M*src$ with $M$
   * being this matrix.
   */
  template <typename block_number>
  void vmult (BlockVector<block_number>       &dst,
              const BlockVector<block_number> &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block column.
   */
  template <typename block_number,
            typename nonblock_number>
  void vmult (BlockVector<block_number>          &dst,
              const Vector<nonblock_number> &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block row.
   */
  template <typename block_number,
            typename nonblock_number>
  void vmult (Vector<nonblock_number>    &dst,
              const BlockVector<block_number> &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block.
   */
  template <typename nonblock_number>
  void vmult (Vector<nonblock_number>       &dst,
              const Vector<nonblock_number> &src) const;

  /**
   * Matrix-vector multiplication:
   * let $dst = M^T*src$ with $M$
   * being this matrix. This
   * function does the same as
   * vmult() but takes the
   * transposed matrix.
   */
  template <typename block_number>
  void Tvmult (BlockVector<block_number>       &dst,
               const BlockVector<block_number> &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block row.
   */
  template <typename block_number,
            typename nonblock_number>
  void Tvmult (BlockVector<block_number>  &dst,
               const Vector<nonblock_number> &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block column.
   */
  template <typename block_number,
            typename nonblock_number>
  void Tvmult (Vector<nonblock_number>    &dst,
               const BlockVector<block_number> &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block.
   */
  template <typename nonblock_number>
  void Tvmult (Vector<nonblock_number>       &dst,
               const Vector<nonblock_number> &src) const;
//@}

  /**
   * @name Preconditioning methods
   */
//@{
  /**
   * Apply the Jacobi
   * preconditioner, which
   * multiplies every element of
   * the <tt>src</tt> vector by the
   * inverse of the respective
   * diagonal element and
   * multiplies the result with the
   * relaxation parameter
   * <tt>omega</tt>.
   *
   * All diagonal blocks must be
   * square matrices for this
   * operation.
   */
  template <class BlockVectorType>
  void precondition_Jacobi (BlockVectorType       &dst,
                            const BlockVectorType &src,
                            const number           omega = 1.) const;

  /**
   * Apply the Jacobi
   * preconditioner to a simple vector.
   *
   * The matrix must be a single
   * square block for this.
   */
  template <typename number2>
  void precondition_Jacobi (Vector<number2>       &dst,
                            const Vector<number2> &src,
                            const number           omega = 1.) const;
//@}

  /**
   * @name Input/Output
   */
//@{
  /**
   * Print the matrix in the usual
   * format, i.e. as a matrix and
   * not as a list of nonzero
   * elements. For better
   * readability, elements not in
   * the matrix are displayed as
   * empty space, while matrix
   * elements which are explicitly
   * set to zero are displayed as
   * such.
   *
   * The parameters allow for a
   * flexible setting of the output
   * format: <tt>precision</tt> and
   * <tt>scientific</tt> are used
   * to determine the number
   * format, where <tt>scientific =
   * false</tt> means fixed point
   * notation.  A zero entry for
   * <tt>width</tt> makes the
   * function compute a width, but
   * it may be changed to a
   * positive value, if output is
   * crude.
   *
   * Additionally, a character for
   * an empty value may be
   * specified.
   *
   * Finally, the whole matrix can
   * be multiplied with a common
   * denominator to produce more
   * readable output, even
   * integers.
   *
   * @attention This function may
   * produce <b>large</b> amounts
   * of output if applied to a
   * large matrix!
   */
  void print_formatted (std::ostream       &out,
                        const unsigned int  precision   = 3,
                        const bool          scientific  = true,
                        const unsigned int  width       = 0,
                        const char         *zero_string = " ",
                        const double        denominator = 1.) const;
//@}
  /** @addtogroup Exceptions
   * @{ */

  /**
   * Exception
   */
  DeclException0 (ExcBlockDimensionMismatch);
  //@}

private:
  /**
   * Pointer to the block sparsity
   * pattern used for this
   * matrix. In order to guarantee
   * that it is not deleted while
   * still in use, we subscribe to
   * it using the SmartPointer
   * class.
   */
  SmartPointer<const BlockSparsityPattern,BlockSparseMatrix<number> > sparsity_pattern;
};



/*@}*/
/* ------------------------- Template functions ---------------------- */



template <typename number>
inline
BlockSparseMatrix<number> &
BlockSparseMatrix<number>::operator = (const double d)
{
  Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());

  for (size_type r=0; r<this->n_block_rows(); ++r)
    for (size_type c=0; c<this->n_block_cols(); ++c)
      this->block(r,c) = d;

  return *this;
}



template <typename number>
template <typename block_number>
inline
void
BlockSparseMatrix<number>::vmult (BlockVector<block_number>       &dst,
                                  const BlockVector<block_number> &src) const
{
  BaseClass::vmult_block_block (dst, src);
}



template <typename number>
template <typename block_number,
          typename nonblock_number>
inline
void
BlockSparseMatrix<number>::vmult (BlockVector<block_number>     &dst,
                                  const Vector<nonblock_number> &src) const
{
  BaseClass::vmult_block_nonblock (dst, src);
}



template <typename number>
template <typename block_number,
          typename nonblock_number>
inline
void
BlockSparseMatrix<number>::vmult (Vector<nonblock_number>         &dst,
                                  const BlockVector<block_number> &src) const
{
  BaseClass::vmult_nonblock_block (dst, src);
}



template <typename number>
template <typename nonblock_number>
inline
void
BlockSparseMatrix<number>::vmult (Vector<nonblock_number>       &dst,
                                  const Vector<nonblock_number> &src) const
{
  BaseClass::vmult_nonblock_nonblock (dst, src);
}



template <typename number>
template <typename block_number>
inline
void
BlockSparseMatrix<number>::Tvmult (BlockVector<block_number>       &dst,
                                   const BlockVector<block_number> &src) const
{
  BaseClass::Tvmult_block_block (dst, src);
}



template <typename number>
template <typename block_number,
          typename nonblock_number>
inline
void
BlockSparseMatrix<number>::Tvmult (BlockVector<block_number>     &dst,
                                   const Vector<nonblock_number> &src) const
{
  BaseClass::Tvmult_block_nonblock (dst, src);
}



template <typename number>
template <typename block_number,
          typename nonblock_number>
inline
void
BlockSparseMatrix<number>::Tvmult (Vector<nonblock_number>         &dst,
                                   const BlockVector<block_number> &src) const
{
  BaseClass::Tvmult_nonblock_block (dst, src);
}



template <typename number>
template <typename nonblock_number>
inline
void
BlockSparseMatrix<number>::Tvmult (Vector<nonblock_number>       &dst,
                                   const Vector<nonblock_number> &src) const
{
  BaseClass::Tvmult_nonblock_nonblock (dst, src);
}



template <typename number>
template <class BlockVectorType>
inline
void
BlockSparseMatrix<number>::
precondition_Jacobi (BlockVectorType       &dst,
                     const BlockVectorType &src,
                     const number           omega) const
{
  Assert (this->n_block_rows() == this->n_block_cols(), ExcNotQuadratic());
  Assert (dst.n_blocks() == this->n_block_rows(),
          ExcDimensionMismatch(dst.n_blocks(), this->n_block_rows()));
  Assert (src.n_blocks() == this->n_block_cols(),
          ExcDimensionMismatch(src.n_blocks(), this->n_block_cols()));

  // do a diagonal preconditioning. uses only
  // the diagonal blocks of the matrix
  for (size_type i=0; i<this->n_block_rows(); ++i)
    this->block(i,i).precondition_Jacobi (dst.block(i),
                                          src.block(i),
                                          omega);
}



template <typename number>
template <typename number2>
inline
void
BlockSparseMatrix<number>::
precondition_Jacobi (Vector<number2>       &dst,
                     const Vector<number2> &src,
                     const number           omega) const
{
  // check number of blocks. the sizes of the
  // single block is checked in the function
  // we call
  Assert (this->n_block_cols() == 1,
          ExcMessage ("This function only works if the matrix has "
                      "a single block"));
  Assert (this->n_block_rows() == 1,
          ExcMessage ("This function only works if the matrix has "
                      "a single block"));

  // do a diagonal preconditioning. uses only
  // the diagonal blocks of the matrix
  this->block(0,0).precondition_Jacobi (dst, src, omega);
}


DEAL_II_NAMESPACE_CLOSE

#endif    // __deal2__block_sparse_matrix_h
