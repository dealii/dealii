//----------------------------  block_sparse_matrix.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_sparse_matrix.h  ---------------------------
#ifndef __deal2__block_sparse_matrix_h
#define __deal2__block_sparse_matrix_h


#include <lac/sparse_matrix.h>
#include <lac/block_sparsity_pattern.h>



template <typename number, int  rows, int columns=rows>
class BlockSparseMatrix : public Subscriptor
{
  public:
				     /**
				      * Type of matrix entries. In analogy to
				      * the STL container classes.
				      */
    typedef number value_type;
    
				     /**
				      * Constructor; initializes the matrix to
				      * be empty, without any structure, i.e.
				      * the matrix is not usable at all. This
				      * constructor is therefore only useful
				      * for matrices which are members of a
				      * class. All other matrices should be
				      * created at a point in the data flow
				      * where all necessary information is
				      * available.
				      *
				      * You have to initialize
				      * the matrix before usage with
				      * #reinit(BlockSparsityPattern)#.
				      */
    BlockSparseMatrix ();

				     /**
				      * Constructor. Takes the given
				      * matrix sparsity structure to
				      * represent the sparsity pattern
				      * of this matrix. You can change
				      * the sparsity pattern later on
				      * by calling the #reinit#
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
				      * long as #reinit# is not called
				      * with a new sparsity structure.
				      */
    BlockSparseMatrix (const BlockSparsityPattern<rows,columns> &sparsity);

    

				     /** 
				      * Pseudo operator only copying
				      * empty objects.
				      */
    BlockSparseMatrix<number,rows,columns> &
    operator = (const BlockSparseMatrix<number,rows,columns> &);


				     /**
				      * Reinitialize the object but
				      * keep to the sparsity pattern
				      * previously used.  This may be
				      * necessary if you #reinit#'d
				      * the sparsity structure and
				      * want to update the size of the
				      * matrix. It only calls #reinit#
				      * on the sub-matrices.
				      */
    virtual void reinit ();

				     /**
				      * Reinitialize the sparse matrix
				      * with the given sparsity
				      * pattern. The latter tells the
				      * matrix how many nonzero
				      * elements there need to be
				      * reserved.
				      *
				      * Basically, this function only
				      * calls #reinit# of the
				      * sub-matrices with the block
				      * sparsity patterns of the
				      * parameter.
				      */
    virtual void reinit (const BlockSparsityPattern<rows,columns> &sparsity);

    
				     /**
				      * Access the block with the
				      * given coordinates.
				      */
    SparseMatrix<number> &
    block (const unsigned int row,
	   const unsigned int column);
    
    
				     /**
				      * Access the block with the
				      * given coordinates. Version for
				      * constant objects.
				      */
    const SparseMatrix<number> &
    block (const unsigned int row,
	   const unsigned int column) const;    

				     /**
				      * Release all memory and return
				      * to a state just like after
				      * having called the default
				      * constructor. It also forgets
				      * the sparsity pattern it was
				      * previously tied to.
				      *
				      * This calls #clear# on all
				      * sub-matrices.
				      */
    virtual void clear ();
    
				     /**
				      * Return whether the object is
				      * empty. It is empty if either
				      * both dimensions are zero or no
				      * #SparsityPattern# is
				      * associated.
				      */
    bool empty () const;

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
				      * Return the number of nonzero
				      * elements of this
				      * matrix. Actually, it returns
				      * the number of entries in the
				      * sparsity pattern; if any of
				      * the entries should happen to
				      * be zero, it is counted anyway.
				      */
    unsigned int n_nonzero_elements () const;
    
				     /**
				      * Set the element #(i,j)# to #value#.
				      * Throws an error if the entry does
				      * not exist. Still, it is allowed to store
				      * zero values in non-existent fields.
				      */
    void set (const unsigned int i,
	      const unsigned int j,
	      const number value);
    
				     /**
				      * Add #value# to the element
				      * #(i,j)#.  Throws an error if
				      * the entry does not
				      * exist. Still, it is allowed to
				      * store zero values in
				      * non-existent fields.
				      */
    void add (const unsigned int i, const unsigned int j,
	      const number value);

				     /**
				      * Copy the given matrix to this
				      * one.  The operation throws an
				      * error if the sparsity patterns
				      * of the two involved matrices
				      * do not point to the same
				      * object, since in this case the
				      * copy operation is
				      * cheaper. Since this operation
				      * is notheless not for free, we
				      * do not make it available
				      * through #operator =#, since
				      * this may lead to unwanted
				      * usage, e.g. in copy arguments
				      * to functions, which should
				      * really be arguments by
				      * reference.
				      *
				      * The source matrix may be a
				      * matrix of arbitrary type, as
				      * long as its data type is
				      * convertible to the data type
				      * of this matrix.
				      *
				      * The function returns a
				      * reference to #this#.
				      */
    template <typename somenumber>
    BlockSparseMatrix<number,rows,columns> &
    copy_from (const BlockSparseMatrix<somenumber,rows,columns> &source);

				     /**
				      * Add #matrix# scaled by
				      * #factor# to this matrix. The
				      * function throws an error if
				      * the sparsity patterns of the
				      * two involved matrices do not
				      * point to the same object,
				      * since in this case the
				      * operation is cheaper.
				      *
				      * The source matrix may be a matrix
				      * of arbitrary type, as long as its
				      * data type is convertible to the
				      * data type of this matrix.
				      */
    template <typename somenumber>
    void add_scaled (const number factor,
		     const BlockSparseMatrix<somenumber,rows,columns> &matrix);
    
				     /**
				      * Return the value of the entry (i,j).
				      * This may be an expensive operation
				      * and you should always take care
				      * where to call this function.
				      * In order to avoid abuse, this function
				      * throws an exception if the wanted
				      * element does not exist in the matrix.
				      */
    number operator () (const unsigned int i, const unsigned int j) const;


				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M*src$ with $M$
				      * being this matrix.
				      */
    template <typename somenumber>
    void vmult (BlockVector<rows,somenumber>          &dst,
		const BlockVector<columns,somenumber> &src) const;
    
				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M^T*src$ with $M$
				      * being this matrix. This
				      * function does the same as
				      * #vmult# but takes the
				      * transposed matrix.
				      */
    template <typename somenumber>
    void Tvmult (BlockVector<columns,somenumber>    &dst,
		 const BlockVector<rows,somenumber> &src) const;
  
				     /**
				      * Adding Matrix-vector
				      * multiplication. Add $M*src$ on
				      * $dst$ with $M$ being this
				      * matrix.
				      */
    template <typename somenumber>
    void vmult_add (BlockVector<rows,somenumber>          &dst,
		    const BlockVector<columns,somenumber> &src) const;
    
				     /**
				      * Adding Matrix-vector
				      * multiplication. Add $M^T*src$
				      * to $dst$ with $M$ being this
				      * matrix. This function does the
				      * same as #vmult_add# but takes
				      * the transposed matrix.
				      */
    template <typename somenumber>
    void Tvmult_add (BlockVector<columns,somenumber>    &dst,
		     const BlockVector<rows,somenumber> &src) const;
  
				     /**
				      * Return the norm of the vector
				      * $v$ with respect to the norm
				      * induced by this matrix,
				      * i.e. $\left(v,Mv\right)$. This
				      * is useful, e.g. in the finite
				      * element context, where the
				      * $L_2$ norm of a function
				      * equals the matrix norm with
				      * respect to the mass matrix of
				      * the vector representing the
				      * nodal values of the finite
				      * element function. Note that
				      * even though the function's
				      * name might suggest something
				      * different, for historic
				      * reasons not the norm but its
				      * square is returned, as defined
				      * above by the scalar product.
				      *
				      * Obviously, the matrix needs to
				      * be square for this operation.
				      */
    template <typename somenumber>
    somenumber
    matrix_norm_square (const BlockVector<rows,somenumber> &v) const;

				     /**
				      * Compute the matrix scalar
				      * product $\left(u,Mv\right)$.
				      */    
    template <typename somenumber>
    somenumber
    matrix_scalar_product (const BlockVector<rows,somenumber>    &u,
			   const BlockVector<columns,somenumber> &v) const;
    
				     /**
				      * Return a (constant) reference
				      * to the underlying sparsity
				      * pattern of this matrix.
				      *
				      * Though the return value is
				      * declared #const#, you should
				      * be aware that it may change if
				      * you call any nonconstant
				      * function of objects which
				      * operate on it.
				      */
    const BlockSparsityPattern<rows,columns> &
    get_sparsity_pattern () const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotBlockSquare);
    
  private:
    				     /**
				      * Pointer to the block sparsity
				      * pattern used for this
				      * matrix. In order to guarantee
				      * that it is not deleted while
				      * still in use, we subscribe to
				      * it using the #SmartPointer#
				      * class.
				      */
    SmartPointer<const BlockSparsityPattern<rows,columns> > sparsity_pattern;

				     /**
				      * Array of sub-matrices.
				      */
    SparseMatrix<number> sub_objects[rows][columns];
};




/* ------------------------- Template functions ---------------------- */


template <typename number, int  rows, int columns>
BlockSparseMatrix<number,rows,columns>::BlockSparseMatrix () :
		sparsity_pattern (0)
{};



template <typename number, int  rows, int columns>
BlockSparseMatrix<number,rows,columns>::
BlockSparseMatrix (const BlockSparsityPattern<rows,columns> &sparsity)
{
  reinit (sparsity);
};



template <typename number, int  rows, int columns>
BlockSparseMatrix<number,rows,columns> &
BlockSparseMatrix<number,rows,columns>::
operator = (const BlockSparseMatrix<number,rows,columns> &m) 
{
  *sparsity_pattern = *m.sparsity_pattern;
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c) = m.block(r,c);
};

 

template <typename number, int  rows, int columns>
void
BlockSparseMatrix<number,rows,columns>::reinit ()
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).reinit ();
};



template <typename number, int  rows, int columns>
void
BlockSparseMatrix<number,rows,columns>::reinit (const BlockSparsityPattern<rows,columns> &sparsity)
{
  sparsity_pattern = &sparsity;
  
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).reinit (sparsity.block(r,c));
};




template <typename number, int rows, int columns>
inline
SparseMatrix<number> &
BlockSparseMatrix<number,rows,columns>::block (const unsigned int row,
					       const unsigned int column)
{
  Assert (row<rows, ExcIndexRange (row, 0, rows));
  Assert (column<columns, ExcIndexRange (column, 0, columns));
  
  return sub_objects[row][column];
};



template <typename number, int rows, int columns>
inline
const SparseMatrix<number> &
BlockSparseMatrix<number,rows,columns>::block (const unsigned int row,
					       const unsigned int column) const
{
  Assert (row<rows, ExcIndexRange (row, 0, rows));
  Assert (column<columns, ExcIndexRange (column, 0, columns));
  
  return sub_objects[row][column];
};



template <typename number, int  rows, int columns>
void
BlockSparseMatrix<number,rows,columns>::clear () 
{
  sparsity_pattern = 0;
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).clear ();
};



template <typename number, int  rows, int columns>
bool
BlockSparseMatrix<number,rows,columns>::empty () const
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      if (block(r,c).empty () == false)
	return false;
  return true;
};



template <typename number, int  rows, int columns>
unsigned int
BlockSparseMatrix<number,rows,columns>::m () const
{
  return sparsity_pattern->n_rows();
};



template <typename number, int  rows, int columns>
unsigned int
BlockSparseMatrix<number,rows,columns>::n () const
{
  return sparsity_pattern->n_cols();
};



template <typename number, int  rows, int columns>
unsigned int
BlockSparseMatrix<number,rows,columns>::n_nonzero_elements () const
{
  return sparsity_pattern->n_nonzero_elements ();
};



template <typename number, int  rows, int columns>
void
BlockSparseMatrix<number,rows,columns>::set (const unsigned int i,
					     const unsigned int j,
					     const number value)
{
  const pair<unsigned int,unsigned int>
    row_index = sparsity_pattern->row_indices.global_to_local (i),
    col_index = sparsity_pattern->column_indices.global_to_local (j);
  block(row_index.first,col_index.first).set (row_index.second,
					      col_index.second,
					      value);
};



template <typename number, int  rows, int columns>
void
BlockSparseMatrix<number,rows,columns>::add (const unsigned int i,
					     const unsigned int j,
					     const number value)
{
  const pair<unsigned int,unsigned int>
    row_index = sparsity_pattern->row_indices.global_to_local (i),
    col_index = sparsity_pattern->column_indices.global_to_local (j);
  block(row_index.first,col_index.first).add (row_index.second,
					      col_index.second,
					      value);
};



template <typename number, int  rows, int columns>
template <typename somenumber>
BlockSparseMatrix<number,rows,columns> &
BlockSparseMatrix<number,rows,columns>::
copy_from (const BlockSparseMatrix<somenumber,rows,columns> &source)
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).copy_from (source.block(r,c));
  
  return *this;
};



template <typename number, int  rows, int columns>
template <typename somenumber>
void
BlockSparseMatrix<number,rows,columns>::
add_scaled (const number factor,
	    const BlockSparseMatrix<somenumber,rows,columns> &matrix)
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).add_scaled (factor, matrix.block(r,c));
};



template <typename number, int  rows, int columns>
number
BlockSparseMatrix<number,rows,columns>::operator () (const unsigned int i,
						     const unsigned int j) const
{
  const pair<unsigned int,unsigned int>
    row_index = sparsity_pattern->row_indices.global_to_local (i),
    col_index = sparsity_pattern->column_indices.global_to_local (j);
  return block(row_index.first,col_index.first) (row_index.second,
						 col_index.second);
};



template <typename number, int  rows, int columns>
template <typename somenumber>
void
BlockSparseMatrix<number,rows,columns>::
vmult (BlockVector<rows,somenumber>          &dst,
       const BlockVector<columns,somenumber> &src) const
{
  for (unsigned int row=0; row<rows; ++row)
    {
      block(row,0).vmult (dst.block(row),
			  src.block(0));
      for (unsigned int col=1; col<columns; ++col)
      block(row,col).vmult_add (dst.block(row),
				src.block(col));
    };
};



template <typename number, int  rows, int columns>
template <typename somenumber>
void
BlockSparseMatrix<number,rows,columns>::
vmult_add (BlockVector<rows,somenumber>          &dst,
	   const BlockVector<columns,somenumber> &src) const
{
  for (unsigned int row=0; row<rows; ++row)
    {
      block(row,0).vmult_add (dst.block(row),
			      src.block(0));
      for (unsigned int col=1; col<columns; ++col)
      block(row,col).vmult_add (dst.block(row),
				src.block(col));
    };
};




template <typename number, int  rows, int columns>
template <typename somenumber>
void
BlockSparseMatrix<number,rows,columns>::
Tvmult (BlockVector<columns,somenumber>   &/*dst*/,
	const BlockVector<rows,somenumber> &/*src*/) const
{
				   // presently not
				   // implemented. should be simple,
				   // but don't have the time right
				   // now.
  Assert (false, ExcNotImplemented());
};



template <typename number, int  rows, int columns>
template <typename somenumber>
void
BlockSparseMatrix<number,rows,columns>::
Tvmult_add (BlockVector<columns,somenumber>    &/*dst*/,
	    const BlockVector<rows,somenumber> &/*src*/) const
{
				   // presently not
				   // implemented. should be simple,
				   // but don't have the time right
				   // now.
  Assert (false, ExcNotImplemented());
};



template <typename number, int  rows, int columns>
template <typename somenumber>
somenumber
BlockSparseMatrix<number,rows,columns>::
matrix_norm_square (const BlockVector<rows,somenumber> &v) const
{
  Assert (rows == columns, ExcMatrixNotBlockSquare());

  somenumber norm_sqr = 0;
  for (unsigned int row=0; row<rows; ++row)
    for (unsigned int col=0; col<columns; ++col)
      if (row==col)
	norm_sqr += block(row,col).matrix_norm_square (v.block(row));
      else
	norm_sqr += block(row,col).matrix_scalar_product (v.block(row),
							  v.block(col));
  return norm_sqr;
};



template <typename number, int  rows, int columns>
template <typename somenumber>
somenumber
BlockSparseMatrix<number,rows,columns>::
matrix_scalar_product (const BlockVector<rows,somenumber>    &u,
		       const BlockVector<columns,somenumber> &v) const
{
  somenumber result = 0;
  for (unsigned int row=0; row<rows; ++row)
    for (unsigned int col=0; col<columns; ++col)
      result += block(row,col).matrix_scalar_product (u.block(row),
						      v.block(col));
  return result;
};



template <typename number, int  rows, int columns>
const BlockSparsityPattern<rows,columns> &
BlockSparseMatrix<number,rows,columns>::get_sparsity_pattern () const
{
  return *sparsity_pattern;
};





#endif    // __deal2__block_sparse_matrix_h
