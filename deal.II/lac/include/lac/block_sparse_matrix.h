//----------------------------  block_sparse_matrix.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004 by the deal.II authors
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


template <typename> class Vector;
template <typename> class BlockVector;


/*! @addtogroup Matrix1
 *@{
 */

/**
 * Blocked sparse matrix. The behaviour of objects of this type is
 * almost as for the SparseMatrix objects, with most of the
 * functions being implemented in both classes. The main difference is
 * that the matrix represented by this object is composed of an array
 * of sparse matrices (i.e. of type SparseMatrix<number>) and all
 * accesses to the elements of this object are relayed to accesses of
 * the base matrices.
 *
 * In addition to the usual matrix access and linear algebra
 * functions, there are functions block() which allow access to the
 * different blocks of the matrix. This may, for example, be of help
 * when you want to implement Schur complement methods, or block
 * preconditioners, where each block belongs to a specific component
 * of the equation you are presently discretizing.
 *
 * Note that the numbers of blocks and rows are implicitly determined
 * by the sparsity pattern objects used.
 *
 * @ref Instantiations: some (<tt>@<float@> @<double@></tt>)
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
	Accessor (const BlockSparseMatrix<number> *m,
		  const unsigned int               row,
		  const unsigned short             index);
	
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

					 /**
					  * Block row of the
					  * element represented by
					  * this object.
					  */
	unsigned int block_row() const;
	
					 /**
					  * Block column of the
					  * element represented by
					  * this object.
					  */
	unsigned int block_column() const;
	
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
	unsigned int row_block;
	
					 /**
					  * First row of block.
					  */
	unsigned int row_start;	
	
					 /**
					  * Number of block column where column lies in.
					  */
	unsigned int col_block;
	
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
					    * Inverse of operator==().
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
				      * necessary if you reinitialized
				      * the sparsity structure and
				      * want to update the size of the
				      * matrix. It only calls
				      * SparseMatrix::reinit() on the
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
				      * calls SparseMatrix::reinit() of the
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
				      * This calls SparseMatrix::clear on all
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
				      * BlockSparsityPattern is
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
				      * Set the element <tt>(i,j)</tt>
				      * to <tt>value</tt>. Throws an
				      * error if the entry does not
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
				      * Add <tt>value</tt> to the element
				      * <tt>(i,j)</tt>.  Throws an error if
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
				      * through operator=(), since
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
				      * reference to <tt>this</tt>.
				      */
    template <typename somenumber>
    BlockSparseMatrix<number> &
    copy_from (const BlockSparseMatrix<somenumber> &source);

				     /**
				      * Add <tt>matrix</tt> scaled by
				      * <tt>factor</tt> to this matrix,
				      * i.e. the matrix <tt>factor*matrix</tt>
				      * is added to <tt>this</tt>. This
				      * function throws an error if the
				      * sparsity patterns of the two involved
				      * matrices do not point to the same
				      * object, since in this case the
				      * operation is cheaper.
				      *
				      * The source matrix may be a sparse
				      * matrix over an arbitrary underlying
				      * scalar type, as long as its data type
				      * is convertible to the data type of
				      * this matrix.
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
				      * operator()() in that it
				      * returns the value of the
				      * matrix entry <tt>(i,j)</tt>. The only
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
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block column.
				      */
    template <typename somenumber>
    void vmult (BlockVector<somenumber>  &dst,
		const Vector<somenumber> &src) const;

				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block row.
				      */
    template <typename somenumber>
    void vmult (Vector<somenumber>            &dst,
		const BlockVector<somenumber> &src) const;

				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block.
				      */
    template <typename somenumber>
    void vmult (Vector<somenumber>       &dst,
		const Vector<somenumber> &src) const;
    
				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M^T*src$ with $M$
				      * being this matrix. This
				      * function does the same as
				      * vmult() but takes the
				      * transposed matrix.
				      */
    template <typename somenumber>
    void Tvmult (BlockVector<somenumber>       &dst,
		 const BlockVector<somenumber> &src) const;
  
				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block row.
				      */
    template <typename somenumber>
    void Tvmult (BlockVector<somenumber>  &dst,
		 const Vector<somenumber> &src) const;

				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block column.
				      */
    template <typename somenumber>
    void Tvmult (Vector<somenumber>            &dst,
		 const BlockVector<somenumber> &src) const;

				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block.
				      */
    template <typename somenumber>
    void Tvmult (Vector<somenumber>       &dst,
		 const Vector<somenumber> &src) const;
    
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
				      * multiplication. Add
				      * <i>M<sup>T</sup>src</i> to
				      * <i>dst</i> with <i>M</i> being
				      * this matrix. This function
				      * does the same as vmult_add()
				      * but takes the transposed
				      * matrix.
				      */
    template <typename somenumber>
    void Tvmult_add (BlockVector<somenumber>       &dst,
		     const BlockVector<somenumber> &src) const;
  
				     /**
				      * Return the norm of the vector
				      * <i>v</i> with respect to the
				      * norm induced by this matrix,
				      * i.e. <i>v<sup>T</sup>Mv)</i>. This
				      * is useful, e.g. in the finite
				      * element context, where the
				      * <i>L<sup>T</sup></i>-norm of a
				      * function equals the matrix
				      * norm with respect to the mass
				      * matrix of the vector
				      * representing the nodal values
				      * of the finite element
				      * function. Note that even
				      * though the function's name
				      * might suggest something
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
				      * Compute the residual
				      * <i>r=b-Ax</i>. Write the
				      * residual into <tt>dst</tt>.
				      */
    template <typename somenumber>
    somenumber residual (BlockVector<somenumber>       &dst,
			 const BlockVector<somenumber> &x,
			 const BlockVector<somenumber> &b) const;

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
    template <typename somenumber>
    void precondition_Jacobi (BlockVector<somenumber>       &dst,
			      const BlockVector<somenumber> &src,
			      const number                   omega = 1.) const;

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
				      *
				      * In contrast to the previous
				      * function, input and output
				      * vectors are not block
				      * vectors. This operation is
				      * thus only allowed if the block
				      * matrix consists of only one
				      * block.
				      */
    template <typename somenumber>
    void precondition_Jacobi (Vector<somenumber>       &dst,
			      const Vector<somenumber> &src,
			      const number             omega = 1.) const;
    
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
				      * first entry of row <tt>r</tt>.
				      */
    const_iterator begin (unsigned int r) const;

				     /**
				      * Final iterator of row <tt>r</tt>.
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
				      * reinit() function.
				      */
    unsigned int rows;

				     /**
				      * Number of block columns in this
				      * matrix. Is set by default to
				      * zero, and is only changed if a
				      * sparsity pattern is given to
				      * the constructor or the
				      * reinit() function.
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



/*@}*/
/* ------------------------- Template functions ---------------------- */

template <typename number>
inline
BlockSparseMatrix<number>::Accessor::
Accessor (const BlockSparseMatrix<number> *matrix,
          const unsigned int               r,
          const unsigned short             i)
		:
                matrix(matrix),
                base_iterator(matrix->block(0,0).begin()),
		row_block(0),
		row_start(0),
		col_block(0),
		col_start(0),
		a_index(0)
{
  Assert (i==0, ExcNotImplemented());

  if (r < matrix->m())
    {
      std::pair<unsigned int,unsigned int> indices
	= matrix->sparsity_pattern->get_row_indices().global_to_local(r);
      row_block = indices.first;
      base_iterator = matrix->block(indices.first, 0).begin(indices.second);
      row_start = matrix->sparsity_pattern
		  ->get_row_indices().local_to_global(row_block, 0);
    }
  else
    {
      row_block = matrix->n_block_rows();
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
unsigned int
BlockSparseMatrix<number>::Accessor::block_row() const
{
  return row_block;
}


template <typename number>
inline
unsigned int
BlockSparseMatrix<number>::Accessor::block_column() const
{
  return col_block;
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
BlockSparseMatrix<number>::const_iterator::
const_iterator(const BlockSparseMatrix<number>* m,
               unsigned int r,
               unsigned short i)
		:
                BlockSparseMatrix<number>::Accessor(m, r, i)
{}



template <typename number>
inline
typename BlockSparseMatrix<number>::const_iterator&
BlockSparseMatrix<number>::const_iterator::operator++ ()
{
  Assert (this->row_block<this->matrix->n_block_rows(), ExcIteratorPastEnd());

				   // Remeber current row inside block
  unsigned int local_row = this->base_iterator->row();
				   // Advance inside block
  ++this->base_iterator;
  ++this->a_index;
				   // If end of row inside block,
				   // advance to next block
  if (this->base_iterator == this->matrix->block(this->row_block, this->col_block).end(local_row))
    {
      if (this->col_block<this->matrix->n_block_cols()-1)
	{
					   // Advance to next block in
					   // row
	  ++this->col_block;
	  this->col_start = this->matrix->sparsity_pattern
		      ->get_column_indices().local_to_global(this->col_block, 0);
	}
      else
	{
					   // Advance to first block
					   // in next row
	  this->col_block = 0;
	  this->col_start = 0;
	  this->a_index = 0;
	  ++local_row;
	  if (local_row>=this->matrix->block(this->row_block,0).m())
	    {
					       // If final row in
					       // block, go to next
					       // block row
	      local_row = 0;
	      ++this->row_block;
	      if (this->row_block < this->matrix->n_block_rows())
		this->row_start = this->matrix->sparsity_pattern
			    ->get_row_indices().local_to_global(this->row_block, 0);
	    }
	}
				       // Finally, set base_iterator
				       // to start of row determined
				       // above
      if (this->row_block < this->matrix->n_block_rows())
	this->base_iterator = this->matrix->block(this->row_block, this->col_block).begin(local_row);
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
BlockSparseMatrix<number>::const_iterator::
operator == (const const_iterator& i) const
{
  if (this->matrix != i->matrix)
    return false;
  
  if (this->row_block == i->row_block
      && this->col_block == i->col_block
      && this->base_iterator == i->base_iterator)
    return true;
  return false;
}



template <typename number>
inline
bool
BlockSparseMatrix<number>::const_iterator::
operator != (const const_iterator& i) const
{
  return !(*this == i);
}



template <typename number>
inline
bool
BlockSparseMatrix<number>::const_iterator::
operator < (const const_iterator& i) const
{
  if (this->row_block<i->row_block)
    return true;
  if (this->row_block == i->row_block)
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
BlockSparseMatrix<number>::vmult (Vector<somenumber>            &dst,
				  const BlockVector<somenumber> &src) const
{
  Assert (rows == 1,
	  ExcDimensionMismatch(1, rows));
  Assert (src.n_blocks() == columns,
	  ExcDimensionMismatch(src.n_blocks(), columns));

  block(0,0).vmult (dst, src.block(0));
  for (unsigned int col=1; col<columns; ++col)
    block(0,col).vmult_add (dst, src.block(col));
}



template <typename number>
template <typename somenumber>
void
BlockSparseMatrix<number>::vmult (BlockVector<somenumber>  &dst,
				  const Vector<somenumber> &src) const
{
  Assert (dst.n_blocks() == rows,
	  ExcDimensionMismatch(dst.n_blocks(), rows));
  Assert (1 == columns,
	  ExcDimensionMismatch(1, columns));

  for (unsigned int row=0; row<rows; ++row)
    block(row,0).vmult (dst.block(row),
			src);
}



template <typename number>
template <typename somenumber>
void
BlockSparseMatrix<number>::vmult (Vector<somenumber>       &dst,
				  const Vector<somenumber> &src) const
{
  Assert (1 == rows,
	  ExcDimensionMismatch(1, rows));
  Assert (1 == columns,
	  ExcDimensionMismatch(1, columns));

  block(0,0).vmult (dst, src);
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
BlockSparseMatrix<number>::Tvmult (BlockVector<somenumber>       &dst,
				   const BlockVector<somenumber> &src) const
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
BlockSparseMatrix<number>::Tvmult (BlockVector<somenumber>  &dst,
				   const Vector<somenumber> &src) const
{
  Assert (dst.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert (1 == n_block_rows(),
	  ExcDimensionMismatch(1, n_block_rows()));

  dst = 0.;
  
  for (unsigned int col=0; col<n_block_cols(); ++col)
    block(0,col).Tvmult_add (dst.block(col), src);
}



template <typename number>
template <typename somenumber>
void
BlockSparseMatrix<number>::Tvmult (Vector<somenumber>            &dst,
				   const BlockVector<somenumber> &src) const
{
  Assert (1 == n_block_cols(),
	  ExcDimensionMismatch(1, n_block_cols()));
  Assert (src.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  block(0,0).Tvmult (dst, src.block(0));
  
  for (unsigned int row=1; row<n_block_rows(); ++row)
    block(row,0).Tvmult_add (dst, src.block(row));
}



template <typename number>
template <typename somenumber>
void
BlockSparseMatrix<number>::Tvmult (Vector<somenumber>       &dst,
				   const Vector<somenumber> &src) const
{
  Assert (1 == n_block_cols(),
	  ExcDimensionMismatch(1, n_block_cols()));
  Assert (1 == n_block_rows(),
	  ExcDimensionMismatch(1, n_block_rows()));

  block(0,0).Tvmult (dst, src);
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
				   // r_i = b_i - \sum_j A_ij x_j.
				   // this can be written as
				   // r_i = b_i - A_i0 x_0 - \sum_{j>0} A_ij x_j.
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
precondition_Jacobi (BlockVector<somenumber>       &dst,
		     const BlockVector<somenumber> &src,
		     const number                  omega) const
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
template <typename somenumber>
void
BlockSparseMatrix<number>::
precondition_Jacobi (Vector<somenumber>       &dst,
		     const Vector<somenumber> &src,
		     const number             omega) const
{
  Assert (rows == columns, ExcMatrixNotBlockSquare());
  Assert (1 == rows,
	  ExcDimensionMismatch(1, rows));
  Assert (1 == columns,
	  ExcDimensionMismatch(1, columns));
  
  block(0,0).precondition_Jacobi (dst, src, omega);
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
