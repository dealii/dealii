//----------------------------  block_sparse_matrix.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_sparse_matrix.h  ---------------------------
#ifndef __deal2__block_sparse_matrix_h
#define __deal2__block_sparse_matrix_h


#include <base/config.h>
#include <base/table.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparsity_pattern.h>
#include <cmath>


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
 * Note that the number of blocks and rows are implicitly determined
 * by the sparsity pattern objects used.
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
template <typename number>
class BlockSparseMatrix : public Subscriptor
{
  public:
				     /**
				      * Type of matrix entries. In analogy to
				      * the STL container classes.
				      */
    typedef number value_type;

    class const_iterator;
    				     /**
				      * Accessor class for iterators
				      */
    class Accessor
    {
      public:
					 /**
					  * Constructor. Since we use
					  * accessors only for read
					  * access, a const matrix
					  * pointer is sufficient.
					  */
	Accessor (const BlockSparseMatrix<number>*,
		  unsigned int row,
		  unsigned short index);
	
					 /**
					  * Row number of the element
					  * represented by this
					  * object.
					  */
	unsigned int row() const;
	
					 /**
					  * Index in row of the element
					  * represented by this
					  * object.
					  */
	unsigned short index() const;
	
					 /**
					  * Column number of the
					  * element represented by
					  * this object.
					  */
	unsigned int column() const;
	
					 /**
					  * Value of this matrix entry.
					  */
	number value() const;
	
      protected:
					 /**
					  * The matrix accessed.
					  */
	const BlockSparseMatrix<number>* matrix;
	
					 /**
					  * Iterator of the underlying matrix class.
					  */
	typename SparseMatrix<number>::const_iterator base_iterator;
	
					 /**
					  * Number of block where row lies in.
					  */
	unsigned int block_row;
	
					 /**
					  * First row of block.
					  */
	unsigned int row_start;	
	
					 /**
					  * Number of block column where column lies in.
					  */
	unsigned int block_col;
	
					 /**
					  * First column of block.
					  */
	unsigned int col_start;
	
					 /**
					  * Index in whole row.
					  */
	unsigned int a_index;

	friend class const_iterator;
    };
    
				     /**
				      * STL conforming iterator.
				      */
    class const_iterator : private Accessor
      {
	public:
					   /**
					    * Constructor.
					    */ 
	const_iterator(const BlockSparseMatrix<number>*,
		       unsigned int row,
		       unsigned short index);
	  
					   /**
					    * Prefix increment.
					    */
	const_iterator& operator++ ();

					   /**
					    * Postfix increment.
					    */
	const_iterator& operator++ (int);

					   /**
					    * Dereferencing operator.
					    */
	const Accessor& operator* () const;

					   /**
					    * Dereferencing operator.
					    */
	const Accessor* operator-> () const;

					   /**
					    * Comparison. True, if
					    * both iterators point to
					    * the same matrix
					    * position.
					    */
	bool operator == (const const_iterator&) const;
					   /**
					    * Inverse of @p{==}.
					    */
	bool operator != (const const_iterator&) const;

					   /**
					    * Comparison
					    * operator. Result is true
					    * if either the first row
					    * number is smaller or if
					    * the row numbers are
					    * equal and the first
					    * index is smaller.
					    */
	bool operator < (const const_iterator&) const;
      };

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
				      * @p{reinit(BlockSparsityPattern)}. The
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
    BlockSparseMatrix (const BlockSparsityPattern &sparsity);

				     /**
				      * Destructor.
				      */
    virtual ~BlockSparseMatrix ();
    
    

				     /** 
				      * Pseudo operator only copying
				      * empty objects. The sizes of
				      * the block matrices need to be
				      * the same.
				      */
    BlockSparseMatrix<number> &
    operator = (const BlockSparseMatrix<number> &);

				     /**
				      * Reinitialize the object but
				      * keep to the sparsity pattern
				      * previously used.  This may be
				      * necessary if you @p{reinit}'d
				      * the sparsity structure and
				      * want to update the size of the
				      * matrix. It only calls
				      * @p{reinit} on the
				      * sub-matrices. The size of this
				      * matrix is unchanged.
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
    virtual void reinit (const BlockSparsityPattern &sparsity);

    
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
				      * Return the number of blocks in
				      * a column. Returns zero if no
				      * sparsity pattern is presently
				      * associated to this matrix.
				      */
    unsigned int n_block_rows () const;
    
				     /**
				      * Return the number of blocks in
				      * a row. Returns zero if no
				      * sparsity pattern is presently
				      * associated to this matrix.
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
				      * Return the number of actually
				      * nonzero elements. Just counts
				      * the number of actually nonzero
				      * elements of all the blocks.
				      */
    unsigned int n_actually_nonzero_elements () const;
    
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
	      const number value);
    
				     /**
				      * Multiply the entire matrix by a
				      * fixed factor.
				      */
    BlockSparseMatrix & operator *= (const number factor);
    
				     /**
				      * Divide the entire matrix by a
				      * fixed factor.
				      */
    BlockSparseMatrix & operator /= (const number factor);
    
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
    BlockSparseMatrix<number> &
    copy_from (const BlockSparseMatrix<somenumber> &source);

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
				      * The source matrix may be a
				      * matrix of arbitrary type, as
				      * long as its data type is
				      * convertible to the data type
				      * of this matrix.
				      */
    template <typename somenumber>
    void add_scaled (const number factor,
		     const BlockSparseMatrix<somenumber> &matrix);
    
				     /**
				      * Return the value of the entry
				      * (i,j).  This may be an
				      * expensive operation and you
				      * should always take care where
				      * to call this function.  In
				      * order to avoid abuse, this
				      * function throws an exception
				      * if the wanted element does not
				      * exist in the matrix.
				      */
    number operator () (const unsigned int i,
			const unsigned int j) const;

				     /**
				      * This function is mostly like
				      * @p{operator()} in that it
				      * returns the value of the
				      * matrix entry @p{(i,j)}. The only
				      * difference is that if this
				      * entry does not exist in the
				      * sparsity pattern, then instead
				      * of raising an exception, zero
				      * is returned. While this may be
				      * convenient in some cases, note
				      * that it is simple to write
				      * algorithms that are slow
				      * compared to an optimal
				      * solution, since the sparsity
				      * of the matrix is not used.
				      */
    number el (const unsigned int i,
	       const unsigned int j) const;

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
    matrix_scalar_product (const BlockVector<somenumber> &u,
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
    somenumber residual (BlockVector<somenumber>       &dst,
			 const BlockVector<somenumber> &x,
			 const BlockVector<somenumber> &b) const;

				     /**
				      * Apply the Jacobi
				      * preconditioner, which
				      * multiplies every element of
				      * the @p{src} vector by the
				      * inverse of the respective
				      * diagonal element and
				      * multiplies the result with the
				      * relaxation parameter @p{omega}.
				      *
				      * All diagonal blocks must be
				      * square matrices for this
				      * operation.
				      */
    template <typename somenumber>
    void precondition_Jacobi (BlockVector<somenumber>       &dst,
			      const BlockVector<somenumber> &src,
			      const number                   omega = 1.) const;

                                     /* Call print functions for 
				      * the SparseMatrix blocks.
				      */
    void print_formatted (std::ostream       &out,
			  const unsigned int  precision   = 3,
			  const bool          scientific  = true,
			  const unsigned int  width       = 0,
			  const char         *zero_string = " ",
			  const double        denominator = 1.) const;

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
    const BlockSparsityPattern &
    get_sparsity_pattern () const;

				     /**
				      * STL-like iterator with the
				      * first entry.
				      */
    const_iterator begin () const;

				     /**
				      * Final iterator.
				      */
    const_iterator end () const;
    
				     /**
				      * STL-like iterator with the
				      * first entry of row @p{r}.
				      */
    const_iterator begin (unsigned int r) const;

				     /**
				      * Final iterator of row @p{r}.
				      */
    const_iterator end (unsigned int r) const;
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotBlockSquare);

				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotInitialized);

    
  private:
				     /**
				      * Number of block rows in this
				      * matrix. Is set by default to
				      * zero, and is only changed if a
				      * sparsity pattern is given to
				      * the constructor or the
				      * @p{reinit} function.
				      */
    unsigned int rows;

				     /**
				      * Number of block columns in this
				      * matrix. Is set by default to
				      * zero, and is only changed if a
				      * sparsity pattern is given to
				      * the constructor or the
				      * @p{reinit} function.
				      */
    unsigned int columns;
    
    				     /**
				      * Pointer to the block sparsity
				      * pattern used for this
				      * matrix. In order to guarantee
				      * that it is not deleted while
				      * still in use, we subscribe to
				      * it using the @p{SmartPointer}
				      * class.
				      */
    SmartPointer<const BlockSparsityPattern> sparsity_pattern;

				     /**
				      * Array of sub-matrices.
				      */
    Table<2,SmartPointer<SparseMatrix<number> > > sub_objects;

    friend class Accessor;
    friend class const_iterator;
};




/* ------------------------- Template functions ---------------------- */

template <typename number>
inline
BlockSparseMatrix<number>::Accessor::Accessor (
  const BlockSparseMatrix<number>* matrix,
  unsigned int r,
  unsigned short i)
		: matrix(matrix),
		  base_iterator(matrix->block(0,0).begin()),
		  block_row(0),
		  row_start(0),
		  block_col(0),
		  col_start(0),
		  a_index(0)
{
  Assert (i==0, ExcNotImplemented());

  if (r < matrix->m())
    {
      std::pair<unsigned int,unsigned int> indices
	= matrix->sparsity_pattern->get_row_indices().global_to_local(r);
      block_row = indices.first;
      base_iterator = matrix->block(indices.first, 0).begin(indices.second);
      row_start = matrix->sparsity_pattern
		  ->get_row_indices().local_to_global(block_row, 0);
    }
  else
    {
      block_row = matrix->n_block_rows();
      base_iterator = matrix->block(0, 0).begin();
    }
}


template <typename number>
inline
unsigned int
BlockSparseMatrix<number>::Accessor::row() const
{
  return row_start + base_iterator->row();
}


template <typename number>
inline
short unsigned int
BlockSparseMatrix<number>::Accessor::index() const
{
  return a_index;
}


template <typename number>
inline
unsigned int
BlockSparseMatrix<number>::Accessor::column() const
{
  return col_start + base_iterator->column();
}


template <typename number>
inline
number
BlockSparseMatrix<number>::Accessor::value () const
{
  return base_iterator->value();
}


//----------------------------------------------------------------------//


template <typename number>
inline
BlockSparseMatrix<number>::const_iterator::const_iterator(
  const BlockSparseMatrix<number>* m,
  unsigned int r,
  unsigned short i)
		: BlockSparseMatrix<number>::Accessor(m, r, i)
{}



template <typename number>
inline
typename BlockSparseMatrix<number>::const_iterator&
BlockSparseMatrix<number>::const_iterator::operator++ ()
{
  Assert (this->block_row<this->matrix->n_block_rows(), ExcIteratorPastEnd());

				   // Remeber current row inside block
  unsigned int local_row = this->base_iterator->row();
				   // Advance inside block
  ++this->base_iterator;
  ++this->a_index;
				   // If end of row inside block,
				   // advance to next block
  if (this->base_iterator == this->matrix->block(this->block_row, this->block_col).end(local_row))
    {
      if (this->block_col<this->matrix->n_block_cols()-1)
	{
					   // Advance to next block in
					   // row
	  ++this->block_col;
	  this->col_start = this->matrix->sparsity_pattern
		      ->get_column_indices().local_to_global(this->block_col, 0);
	}
      else
	{
					   // Advance to first block
					   // in next row
	  this->block_col = 0;
	  this->col_start = 0;
	  this->a_index = 0;
	  ++local_row;
	  if (local_row>=this->matrix->block(this->block_row,0).m())
	    {
					       // If final row in
					       // block, go to next
					       // block row
	      local_row = 0;
	      ++this->block_row;
	      if (this->block_row < this->matrix->n_block_rows())
		this->row_start = this->matrix->sparsity_pattern
			    ->get_row_indices().local_to_global(this->block_row, 0);
	    }
	}
				       // Finally, set base_iterator
				       // to start of row determined
				       // above
      if (this->block_row < this->matrix->n_block_rows())
	this->base_iterator = this->matrix->block(this->block_row, this->block_col).begin(local_row);
      else
					 // Set base_iterator to a
					 // defined state for
					 // comparison. This is the
					 // end() state.
	this->base_iterator = this->matrix->block(0, 0).begin();
    }
  return *this;
}



//  template <typename number>
//  inline
//  const_iterator&
//  BlockSparseMatrix<number>::const_iterator::operator++ (int)
//  {
//    Assert (false, ExcNotImplemented());
//  }




template <typename number>
inline
const typename BlockSparseMatrix<number>::Accessor&
BlockSparseMatrix<number>::const_iterator::operator* () const
{
  return *this;
}


template <typename number>
inline
const typename BlockSparseMatrix<number>::Accessor*
BlockSparseMatrix<number>::const_iterator::operator-> () const
{
  return this;
}



template <typename number>
inline
bool
BlockSparseMatrix<number>::const_iterator::operator == (const const_iterator& i) const
{
  if (this->matrix != i->matrix)
    return false;
  
  if (this->block_row == i->block_row
      && this->block_col == i->block_col
      && this->base_iterator == i->base_iterator)
    return true;
  return false;
}



template <typename number>
inline
bool
BlockSparseMatrix<number>::const_iterator::operator != (const const_iterator& i) const
{
  return !(*this == i);
}



template <typename number>
inline
bool
BlockSparseMatrix<number>::const_iterator::operator < (const const_iterator& i) const
{
  if (this->block_row<i->block_row)
    return true;
  if (this->block_row == i->block_row)
    {
      if (this->base_iterator->row() < i->base_iterator->row())
	return true;
      if (this->base_iterator->row() == i->base_iterator->row())
	{
	  if (this->a_index < i->a_index)
	    return true;
	}
    }
  return false;
}

//----------------------------------------------------------------------//


template <typename number>
inline
unsigned int
BlockSparseMatrix<number>::n_block_cols () const
{
  return columns;
}



template <typename number>
inline
unsigned int
BlockSparseMatrix<number>::n_block_rows () const
{
  return rows;
}


template <typename number>
inline
SparseMatrix<number> &
BlockSparseMatrix<number>::block (const unsigned int row,
				  const unsigned int column)
{
  Assert (row<rows, ExcIndexRange (row, 0, rows));
  Assert (column<columns, ExcIndexRange (column, 0, columns));
  
  return *sub_objects[row][column];
}



template <typename number>
inline
const SparseMatrix<number> &
BlockSparseMatrix<number>::block (const unsigned int row,
				  const unsigned int column) const
{
  Assert (row<rows, ExcIndexRange (row, 0, rows));
  Assert (column<columns, ExcIndexRange (column, 0, columns));
  
  return *sub_objects[row][column];
}



template <typename number>
inline
unsigned int
BlockSparseMatrix<number>::m () const
{
  return sparsity_pattern->n_rows();
}



template <typename number>
inline
unsigned int
BlockSparseMatrix<number>::n () const
{
  return sparsity_pattern->n_cols();
}



template <typename number>
inline
void
BlockSparseMatrix<number>::set (const unsigned int i,
				const unsigned int j,
				const number value)
{
  const std::pair<unsigned int,unsigned int>
    row_index = sparsity_pattern->row_indices.global_to_local (i),
    col_index = sparsity_pattern->column_indices.global_to_local (j);
  block(row_index.first,col_index.first).set (row_index.second,
					      col_index.second,
					      value);
}



template <typename number>
inline
BlockSparseMatrix<number> &
BlockSparseMatrix<number>::operator *= (const number factor)
{
  Assert (columns != 0, ExcMatrixNotInitialized());
  Assert (rows != 0, ExcMatrixNotInitialized());

  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c) *= factor;

  return *this;
}



template <typename number>
inline
BlockSparseMatrix<number> &
BlockSparseMatrix<number>::operator /= (const number factor)
{
  Assert (columns != 0, ExcMatrixNotInitialized());
  Assert (rows != 0, ExcMatrixNotInitialized());
  Assert (factor !=0, ExcDivideByZero());

  const number factor_inv = 1. / factor;

  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c) *= factor_inv;

  return *this;
}



template <typename number>
inline
void
BlockSparseMatrix<number>::add (const unsigned int i,
				const unsigned int j,
				const number value)
{
  const std::pair<unsigned int,unsigned int>
    row_index = sparsity_pattern->row_indices.global_to_local (i),
    col_index = sparsity_pattern->column_indices.global_to_local (j);
  block(row_index.first,col_index.first).add (row_index.second,
					      col_index.second,
					      value);
}



template <typename number>
inline
number
BlockSparseMatrix<number>::operator () (const unsigned int i,
					const unsigned int j) const
{
  const std::pair<unsigned int,unsigned int>
    row_index = sparsity_pattern->row_indices.global_to_local (i),
    col_index = sparsity_pattern->column_indices.global_to_local (j);
  return block(row_index.first,col_index.first) (row_index.second,
						 col_index.second);
}



template <typename number>
inline
number
BlockSparseMatrix<number>::el (const unsigned int i,
			       const unsigned int j) const
{
  const std::pair<unsigned int,unsigned int>
    row_index = sparsity_pattern->row_indices.global_to_local (i),
    col_index = sparsity_pattern->column_indices.global_to_local (j);
  return block(row_index.first,col_index.first).el (row_index.second,
						    col_index.second);
}




template <typename number>
template <typename somenumber>
BlockSparseMatrix<number> &
BlockSparseMatrix<number>::
copy_from (const BlockSparseMatrix<somenumber> &source)
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).copy_from (source.block(r,c));
  
  return *this;
}



template <typename number>
template <typename somenumber>
void
BlockSparseMatrix<number>::
add_scaled (const number factor,
	    const BlockSparseMatrix<somenumber> &matrix)
{
  for (unsigned int r=0; r<rows; ++r)
    for (unsigned int c=0; c<columns; ++c)
      block(r,c).add_scaled (factor, matrix.block(r,c));
}




template <typename number>
template <typename somenumber>
void
BlockSparseMatrix<number>::vmult (BlockVector<somenumber>       &dst,
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
}



template <typename number>
template <typename somenumber>
void
BlockSparseMatrix<number>::vmult_add (BlockVector<somenumber>       &dst,
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
}




template <typename number>
template <typename somenumber>
void
BlockSparseMatrix<number>::Tvmult (BlockVector<somenumber>& dst,
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
}



template <typename number>
template <typename somenumber>
void
BlockSparseMatrix<number>::Tvmult_add (BlockVector<somenumber>& dst,
				       const BlockVector<somenumber>& src) const
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
}



template <typename number>
template <typename somenumber>
somenumber
BlockSparseMatrix<number>::matrix_norm_square (const BlockVector<somenumber> &v) const
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
}



template <typename number>
template <typename somenumber>
somenumber
BlockSparseMatrix<number>::
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
}



template <typename number>
template <typename somenumber>
somenumber
BlockSparseMatrix<number>::
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
  return std::sqrt(res);
}



template <typename number>
template <typename somenumber>
void
BlockSparseMatrix<number>::
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
}


template <typename number>
inline
typename BlockSparseMatrix<number>::const_iterator
BlockSparseMatrix<number>::begin () const
{
  return const_iterator(this, 0, 0);
}

template <typename number>
inline
typename BlockSparseMatrix<number>::const_iterator
BlockSparseMatrix<number>::end () const
{
  return const_iterator(this, m(), 0);
}

template <typename number>
inline
typename BlockSparseMatrix<number>::const_iterator
BlockSparseMatrix<number>::begin (unsigned int r) const
{
  Assert (r<m(), ExcIndexRange(r,0,m()));
  return const_iterator(this, r, 0);
}

template <typename number>
inline
typename BlockSparseMatrix<number>::const_iterator
BlockSparseMatrix<number>::end (unsigned int r) const
{
  Assert (r<m(), ExcIndexRange(r,0,m()));
  return const_iterator(this, r+1, 0);
}

#endif    // __deal2__block_sparse_matrix_h
