//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__block_sparse_matrix_ez_h
#define __deal2__block_sparse_matrix_ez_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/table.h>
#include <base/smartpointer.h>
#include <lac/block_indices.h>
#include <lac/sparse_matrix_ez.h>

template <typename Number> class BlockVector;

/**
 * A block matrix consisting of blocks of type @p{SparseMatrixEZ}.
 *
 * Like the other Block-objects, this matrix can be used like a
 * @p{SparseMatrixEZ}, when it comes to access to entries. Then, there
 * are functions for the multiplication with @p{BlockVector} and
 * access to the individual blocks.
 *
 * @author Guido Kanschat, 2002, 2003
 */
template<typename Number>
class BlockSparseMatrixEZ : public Subscriptor
{
  public:
				     /**
				      * Default constructor. The
				      * result is an empty object with
				      * zero dimensions.
				      */
    BlockSparseMatrixEZ ();

				     /**
				      * Constructor setting up an
				      * object with given unmber of
				      * block rows and columns. The
				      * blocks themselves still have
				      * zero dimension.
				      */
    BlockSparseMatrixEZ (const unsigned int block_rows,
			 const unsigned int block_cols);

				     /**
				      * Copy constructor. This is
				      * needed for some container
				      * classes. It creates an object
				      * of the same number of block
				      * rows and columns. Since it
				      * calls the copy constructor of
				      * @ref{SparseMatrixEZ}, the
				      * block s must be empty.
				      */
    BlockSparseMatrixEZ (const BlockSparseMatrixEZ<Number>&);

				     /**
				      * Copy operator. Like the copy
				      * constructor, this may be
				      * called for objects with empty
				      * blocks only.
				      */
    BlockSparseMatrixEZ & operator = (const BlockSparseMatrixEZ<Number>&);

				     /**
				      * Set matrix to zero dimensions
				      * and release memory.
				      */
    void clear ();
    
				     /**
				      * Initialize to given block
				      * numbers.  After this
				      * operation, the matrix will
				      * have the block dimensions
				      * provided. Each block will have
				      * zero dimensions and must be
				      * initialized
				      * subsequently. After setting
				      * the sizes of the blocks,
				      * @ref{collect_sizes} must be
				      * called to update internal data
				      * structures.
				      */
    void reinit (const unsigned int n_block_rows,
		 const unsigned int n_block_cols);
				     /**
				      * This function collects the
				      * sizes of the sub-objects and
				      * stores them in internal
				      * arrays, in order to be able to
				      * relay global indices into the
				      * matrix to indices into the
				      * subobjects. You *must* call
				      * this function each time after
				      * you have changed the size of
				      * the sub-objects.
				      */
    void collect_sizes ();
    
				     /**
				      * Access the block with the
				      * given coordinates.
				      */
    SparseMatrixEZ<Number>&
    block (const unsigned int row,
	   const unsigned int column);
    
    
				     /**
				      * Access the block with the
				      * given coordinates. Version for
				      * constant objects.
				      */
    const SparseMatrixEZ<Number>&
    block (const unsigned int row,
	   const unsigned int column) const;    

    				     /**
				      * Return the number of blocks in a
				      * column.
				      */
    unsigned int n_block_rows () const;
    
				     /**
				      * Return the number of blocks in a
				      * row.
				      */
    unsigned int n_block_cols () const;
  
				     /**
				      * Return whether the object is
				      * empty. It is empty if no
				      * memory is allocated, which is
				      * the same as that both
				      * dimensions are zero. This
				      * function is just the
				      * concatenation of the
				      * respective call to all
				      * sub-matrices.
				      */
    bool empty () const;

    				     /**
				      * Return number of rows of this
				      * matrix, which equals the
				      * dimension of the image
				      * space. It is the sum of rows
				      * of the rows of sub-matrices.
				      */
    unsigned int n_rows () const;

				     /**
				      * Return number of columns of
				      * this matrix, which equals the
				      * dimension of the range
				      * space. It is the sum of
				      * columns of the columns of
				      * sub-matrices.
				      */
    unsigned int n_cols () const;

				     /**
				      * Return the dimension of the
				      * image space.  To remember: the
				      * matrix is of dimension
				      * $m \times n$.
				      */
    unsigned int m () const;
    
				     /**
				      * Return the dimension of the
				      * range space.  To remember: the
				      * matrix is of dimension
				      * $m \times n$.
				      */
    unsigned int n () const;

				     /**
				      * Set the element @p{(i,j)} to
				      * @p{value}.  Throws an error if
				      * the entry does not
				      * exist. Still, it is allowed to
				      * store zero values in
				      * non-existent fields.
				      */
    void set (const unsigned int i,
	      const unsigned int j,
	      const Number value);
    
				     /**
				      * Add @p{value} to the element
				      * @p{(i,j)}.  Throws an error if
				      * the entry does not
				      * exist. Still, it is allowed to
				      * store zero values in
				      * non-existent fields.
				      */
    void add (const unsigned int i, const unsigned int j,
	      const Number value);

    
				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M*src$ with $M$
				      * being this matrix.
				      */
    template <typename somenumber>
    void vmult (BlockVector<somenumber>       &dst,
		const BlockVector<somenumber> &src) const;
    
				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M^T*src$ with $M$
				      * being this matrix. This
				      * function does the same as
				      * @p{vmult} but takes the
				      * transposed matrix.
				      */
    template <typename somenumber>
    void Tvmult (BlockVector<somenumber>       &dst,
		 const BlockVector<somenumber> &src) const;
  
				     /**
				      * Adding Matrix-vector
				      * multiplication. Add $M*src$ on
				      * $dst$ with $M$ being this
				      * matrix.
				      */
    template <typename somenumber>
    void vmult_add (BlockVector<somenumber>       &dst,
		    const BlockVector<somenumber> &src) const;
    
				     /**
				      * Adding Matrix-vector
				      * multiplication. Add $M^T*src$
				      * to $dst$ with $M$ being this
				      * matrix. This function does the
				      * same as @p{vmult_add} but takes
				      * the transposed matrix.
				      */
    template <typename somenumber>
    void Tvmult_add (BlockVector<somenumber>       &dst,
		     const BlockVector<somenumber> &src) const;

    
  private:
				     /**
				      * Object storing and managing
				      * the transformation of row
				      * indices to indices of the
				      * sub-objects.
				      */
    BlockIndices    row_indices;
    
				     /**
				      * Object storing and managing
				      * the transformation of column
				      * indices to indices of the
				      * sub-objects.
				      */
    BlockIndices    column_indices;

				     /**
				      * The actual matrices
				      */
    Table<2, SparseMatrixEZ<Number> > blocks;
};

  
/*----------------------------------------------------------------------*/


template <typename Number>
inline
unsigned int
BlockSparseMatrixEZ<Number>::n_block_rows () const
{
  return row_indices.size();
}



template <typename Number>
inline
unsigned int
BlockSparseMatrixEZ<Number>::n_rows () const
{
  return row_indices.total_size();
}



template <typename Number>
inline
unsigned int
BlockSparseMatrixEZ<Number>::n_block_cols () const
{
  return column_indices.size();
}



template <typename Number>
inline
unsigned int
BlockSparseMatrixEZ<Number>::n_cols () const
{
  return column_indices.total_size();
}



template <typename Number>
inline
SparseMatrixEZ<Number> &
BlockSparseMatrixEZ<Number>::block (const unsigned int row,
				    const unsigned int column)
{
  Assert (row<n_block_rows(), ExcIndexRange (row, 0, n_block_rows()));
  Assert (column<n_block_cols(), ExcIndexRange (column, 0, n_block_cols()));
  
  return blocks[row][column];
};



template <typename Number>
inline
const SparseMatrixEZ<Number> &
BlockSparseMatrixEZ<Number>::block (const unsigned int row,
				    const unsigned int column) const
{
  Assert (row<n_block_rows(), ExcIndexRange (row, 0, n_block_rows()));
  Assert (column<n_block_cols(), ExcIndexRange (column, 0, n_block_cols()));
  
  return blocks[row][column];
};



template <typename Number>
inline
unsigned int
BlockSparseMatrixEZ<Number>::m () const
{
  return n_rows();
};



template <typename Number>
inline
unsigned int
BlockSparseMatrixEZ<Number>::n () const
{
  return n_cols();
};



template <typename Number>
inline
void
BlockSparseMatrixEZ<Number>::set (const unsigned int i,
				  const unsigned int j,
				  const Number value)
{
  const std::pair<unsigned int,unsigned int>
    row_index = row_indices.global_to_local (i),
    col_index = column_indices.global_to_local (j);
  block(row_index.first,col_index.first).set (row_index.second,
					      col_index.second,
					      value);
};



template <typename Number>
inline
void
BlockSparseMatrixEZ<Number>::add (const unsigned int i,
				  const unsigned int j,
				  const Number value)
{
  const std::pair<unsigned int,unsigned int>
    row_index = row_indices.global_to_local (i),
    col_index = column_indices.global_to_local (j);
  block(row_index.first,col_index.first).add (row_index.second,
					      col_index.second,
					      value);
};


template <typename Number>
template <typename somenumber>
void
BlockSparseMatrixEZ<Number>::vmult (BlockVector<somenumber>       &dst,
				    const BlockVector<somenumber> &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  dst = 0.;
  
  for (unsigned int row=0; row<n_block_rows(); ++row)
    {
      for (unsigned int col=0; col<n_block_cols(); ++col)
	block(row,col).vmult_add (dst.block(row),
				  src.block(col));
    };
};



template <typename Number>
template <typename somenumber>
void
BlockSparseMatrixEZ<Number>::vmult_add (
  BlockVector<somenumber>       &dst,
  const BlockVector<somenumber> &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  for (unsigned int row=0; row<n_block_rows(); ++row)
    {
      for (unsigned int col=0; col<n_block_cols(); ++col)
	block(row,col).vmult_add (dst.block(row),
				  src.block(col));
    };
};




template <typename Number>
template <typename somenumber>
void
BlockSparseMatrixEZ<Number>::Tvmult (
  BlockVector<somenumber>   &dst,
  const BlockVector<somenumber>& src) const
{
  Assert (dst.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert (src.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  dst = 0.;
  
  for (unsigned int row=0; row<n_block_rows(); ++row)
    {
      for (unsigned int col=0; col<n_block_cols(); ++col)
	block(row,col).Tvmult_add (dst.block(col),
				   src.block(row));
    };
};



template <typename Number>
template <typename somenumber>
void
BlockSparseMatrixEZ<Number>::Tvmult_add (
  BlockVector<somenumber>    &dst,
  const BlockVector<somenumber> &src) const
{
  Assert (dst.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert (src.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  for (unsigned int row=0; row<n_block_rows(); ++row)
    {
      for (unsigned int col=0; col<n_block_cols(); ++col)
	block(row,col).Tvmult_add (dst.block(col),
				   src.block(row));
    };
};


#endif //__deal2__block_sparse_matrix_ez_h
