//----------------------------  block_sparse_matrix.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
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

template <typename Number> class BlockVector;

/**
 * Blocked sparse matrix. The behaviour of objects of this type is
 * almost as for the @p{SparseMatrix<...>} objects, with most of the
 * functions being implemented in both classes. The main difference is
 * that the matrix represented by this object is composed of an array
 * of sparse matrices (i.e. of type @p{SparseMatrix<number>}) and all
 * accesses to the elements of this object are relayed to accesses of
 * the base matrices.
 *
 * In addition to the usual matrix access and linear algebra
 * functions, there are functions @p{block} which allow access to the
 * different blocks of the matrix. This may, for example, be of help
 * when you want to implement Schur complement methods, or block
 * preconditioners, where each block belongs to a specific component
 * of the equation you are presently discretizing.
 *
 *
 * @sect2{On template instantiations}
 *
 * Member functions of this class are either implemented in this file
 * or in a file of the same name with suffix ``.templates.h''. For the
 * most common combinations of the template parameters, instantiations
 * of this class are provided in a file with suffix ``.cc'' in the
 * ``source'' directory. If you need an instantiation that is not
 * listed there, you have to include this file along with the
 * corresponding ``.templates.h'' file and instantiate the respective
 * class yourself.
 *
 * @author Wolfgang Bangerth, 2000
 */
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
				      * @p{reinit(BlockSparsityPattern)}.
				      */
    BlockSparseMatrix ();

				     /**
				      * Constructor. Takes the given
				      * matrix sparsity structure to
				      * represent the sparsity pattern
				      * of this matrix. You can change
				      * the sparsity pattern later on
				      * by calling the @p{reinit}
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
				      * long as @p{reinit} is not called
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
				      * necessary if you @p{reinit}'d
				      * the sparsity structure and
				      * want to update the size of the
				      * matrix. It only calls @p{reinit}
				      * on the sub-matrices.
				      *
				      * If the sparsity pattern has
				      * not changed, then the effect
				      * of this function is simply to
				      * reset all matrix entries to
				      * zero.
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
				      * calls @p{reinit} of the
				      * sub-matrices with the block
				      * sparsity patterns of the
				      * parameter.
				      *
				      * The elements of the matrix are
				      * set to zero by this function.
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
				      * This calls @p{clear} on all
				      * sub-matrices.
				      */
    virtual void clear ();
    
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
				      * empty. It is empty if either
				      * both dimensions are zero or no
				      * @p{SparsityPattern} is
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
				      * Set the element @p{(i,j)} to @p{value}.
				      * Throws an error if the entry does
				      * not exist. Still, it is allowed to store
				      * zero values in non-existent fields.
				      */
    void set (const unsigned int i,
	      const unsigned int j,
	      const number value);
    
				     /**
				      * Add @p{value} to the element
				      * @p{(i,j)}.  Throws an error if
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
				      * through @p{operator =}, since
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
				      * reference to @p{this}.
				      */
    template <typename somenumber>
    BlockSparseMatrix<number,rows,columns> &
    copy_from (const BlockSparseMatrix<somenumber,rows,columns> &source);

				     /**
				      * Add @p{matrix} scaled by
				      * @p{factor} to this matrix. The
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
    void vmult (BlockVector<somenumber>          &dst,
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
    void Tvmult (BlockVector<somenumber>    &dst,
		 const BlockVector<somenumber> &src) const;
  
				     /**
				      * Adding Matrix-vector
				      * multiplication. Add $M*src$ on
				      * $dst$ with $M$ being this
				      * matrix.
				      */
    template <typename somenumber>
    void vmult_add (BlockVector<somenumber>          &dst,
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
    void Tvmult_add (BlockVector<somenumber>    &dst,
		     const BlockVector<somenumber> &src) const;
  
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
    matrix_norm_square (const BlockVector<somenumber> &v) const;

				     /**
				      * Compute the matrix scalar
				      * product $\left(u,Mv\right)$.
				      */    
    template <typename somenumber>
    somenumber
    matrix_scalar_product (const BlockVector<somenumber>    &u,
			   const BlockVector<somenumber> &v) const;
    
				     /**
				      * Compute the residual of an
				      * equation @p{Ax=b}, where the
				      * residual is defined to be
				      * @p{r=b-Ax} with @p{x} typically
				      * being an approximate of the
				      * true solution of the
				      * equation. Write the residual
				      * into @p{dst}.
				      */
    template <typename somenumber>
    somenumber residual (BlockVector<somenumber>          &dst,
			 const BlockVector<somenumber> &x,
			 const BlockVector<somenumber>    &b) const;

				     /**
				      * Apply the Jacobi
				      * preconditioner, which
				      * multiplies every element of
				      * the @p{src} vector by the
				      * inverse of the respective
				      * diagonal element and
				      * multiplies the result with the
				      * damping factor @p{omega}.
				      *
				      * The matrix needs to be square
				      * for this operation.
				      */
    template <typename somenumber>
    void precondition_Jacobi (BlockVector<somenumber>          &dst,
			      const BlockVector<somenumber> &src,
			      const number                           omega = 1.) const;

				     /**
				      * Return a (constant) reference
				      * to the underlying sparsity
				      * pattern of this matrix.
				      *
				      * Though the return value is
				      * declared @p{const}, you should
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
				      * it using the @p{SmartPointer}
				      * class.
				      */
    SmartPointer<const BlockSparsityPattern<rows,columns> > sparsity_pattern;

				     /**
				      * Array of sub-matrices.
				      */
    SparseMatrix<number> sub_objects[rows][columns];
};




/* ------------------------- Template functions ---------------------- */


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
inline
unsigned int
BlockSparseMatrix<number,rows,columns>::m () const
{
  return sparsity_pattern->n_rows();
};



template <typename number, int  rows, int columns>
inline
unsigned int
BlockSparseMatrix<number,rows,columns>::n () const
{
  return sparsity_pattern->n_cols();
};



template <typename number, int  rows, int columns>
inline
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
inline
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
inline
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
template <typename somenumber>
void
BlockSparseMatrix<number,rows,columns>::
vmult (BlockVector<somenumber>          &dst,
       const BlockVector<somenumber> &src) const
{
  Assert (dst.n_blocks() == rows,
	  ExcDimensionMismatch(dst.n_blocks(), rows));
  Assert (src.n_blocks() == columns,
	  ExcDimensionMismatch(src.n_blocks(), columns));

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
vmult_add (BlockVector<somenumber>          &dst,
	   const BlockVector<somenumber> &src) const
{
  Assert (dst.n_blocks() == rows,
	  ExcDimensionMismatch(dst.n_blocks(), rows));
  Assert (src.n_blocks() == columns,
	  ExcDimensionMismatch(src.n_blocks(), columns));

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
Tvmult (BlockVector<somenumber>   &/*dst*/,
	const BlockVector<somenumber> &/*src*/) const
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
Tvmult_add (BlockVector<somenumber>    &/*dst*/,
	    const BlockVector<somenumber> &/*src*/) const
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
matrix_norm_square (const BlockVector<somenumber> &v) const
{
  Assert (rows == columns, ExcMatrixNotBlockSquare());
  Assert (v.n_blocks() == rows,
	  ExcDimensionMismatch(v.n_blocks(), rows));

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
matrix_scalar_product (const BlockVector<somenumber>    &u,
		       const BlockVector<somenumber> &v) const
{
  Assert (u.n_blocks() == rows,
	  ExcDimensionMismatch(u.n_blocks(), rows));
  Assert (v.n_blocks() == columns,
	  ExcDimensionMismatch(v.n_blocks(), columns));

  somenumber result = 0;
  for (unsigned int row=0; row<rows; ++row)
    for (unsigned int col=0; col<columns; ++col)
      result += block(row,col).matrix_scalar_product (u.block(row),
						      v.block(col));
  return result;
};



template <typename number, int  rows, int columns>
template <typename somenumber>
somenumber
BlockSparseMatrix<number,rows,columns>::
residual (BlockVector<somenumber>          &dst,
	  const BlockVector<somenumber> &x,
	  const BlockVector<somenumber>    &b) const
{
  Assert (dst.n_blocks() == rows,
	  ExcDimensionMismatch(dst.n_blocks(), rows));
  Assert (b.n_blocks() == rows,
	  ExcDimensionMismatch(b.n_blocks(), rows));
  Assert (x.n_blocks() == columns,
	  ExcDimensionMismatch(x.n_blocks(), columns));
				   // in block notation, the residual is
				   // @p{r_i = b_i - \sum_j A_ij x_j}.
				   // this can be written as
				   // @p{r_i = b_i - A_i0 x_0 - \sum_{j>0} A_ij x_j}.
				   //
				   // for the first two terms, we can
				   // call the residual function of
				   // A_i0. for the other terms, we
				   // use vmult_add. however, we want
				   // to subtract, so in order to
				   // avoid a temporary vector, we
				   // perform a sign change of the
				   // first two term before, and after
				   // adding up
  for (unsigned int row=0; row<rows; ++row)
    {
      block(row,0).residual (dst.block(row),
			     x.block(0),
			     b.block(row));
      
      for (unsigned int i=0; i<dst.block(row).size(); ++i)
	dst.block(row)(i) = -dst.block(row)(i);
      
      for (unsigned int col=1; col<columns; ++col)
	block(row,col).vmult_add (dst.block(row),
				  x.block(col));

      for (unsigned int i=0; i<dst.block(row).size(); ++i)
	dst.block(row)(i) = -dst.block(row)(i);
    };

  somenumber res = 0;
  for (unsigned int row=0; row<rows; ++row)
    res += dst.block(row).norm_sqr ();
  return sqrt(res);
};



template <typename number, int  rows, int columns>
template <typename somenumber>
void
BlockSparseMatrix<number,rows,columns>::
precondition_Jacobi (BlockVector<somenumber>          &dst,
		     const BlockVector<somenumber> &src,
		     const number                           omega) const
{
  Assert (rows == columns, ExcMatrixNotBlockSquare());
  Assert (dst.n_blocks() == rows,
	  ExcDimensionMismatch(dst.n_blocks(), rows));
  Assert (src.n_blocks() == columns,
	  ExcDimensionMismatch(src.n_blocks(), columns));
  
  for (unsigned int i=0; i<rows; ++i)
    block(i,i).precondition_Jacobi (dst.block(i),
				    src.block(i),
				    omega);
};


template <typename number, int rows, int columns>
inline
unsigned int
BlockSparseMatrix<number,rows,columns>::n_block_cols () const
{
  return columns;
}



template <typename number, int rows, int columns>
inline
unsigned int
BlockSparseMatrix<number,rows,columns>::n_block_rows () const
{
  return rows;
}


#endif    // __deal2__block_sparse_matrix_h
