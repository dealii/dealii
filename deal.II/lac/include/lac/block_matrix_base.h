//----------------------------  block_matrix_base.h  ---------------------------
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
//----------------------------  block_matrix_base.h  ---------------------------
#ifndef __deal2__block_matrix_base_h
#define __deal2__block_matrix_base_h


#include <base/config.h>
#include <base/table.h>
#include <base/smartpointer.h>
#include <lac/block_indices.h>

#include <cmath>



/*! @addtogroup Matrix1
 *@{
 */


namespace internal
{

/**
 * Namespace in which iterators in block matrices are implemented.
 *
 * @author Wolfgang Bangerth, 2004
 */
  namespace BlockMatrixIterators
  {
    template <typename> class ConstIterator;
    
    				     /**
				      * Accessor class for iterators
				      */
    template <class BlockMatrix>
    class Accessor
    {
      public:
                                         /**
                                          * Typedef the value type of the
                                          * matrix we point into.
                                          */
        typedef typename BlockMatrix::value_type value_type;
        
					 /**
					  * Constructor. Since we use
					  * accessors only for read
					  * access, a const matrix
					  * pointer is sufficient.
					  */
	Accessor (const BlockMatrix  *m,
		  const unsigned int  row,
		  const unsigned int  index);
	
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
	unsigned int index() const;
	
					 /**
					  * Column number of the
					  * element represented by
					  * this object.
					  */
	unsigned int column() const;
	
					 /**
					  * Value of this matrix entry.
					  */
	value_type value() const;

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
	const BlockMatrix * matrix;
	
					 /**
					  * Iterator of the underlying matrix
					  * class.
					  */
	typename BlockMatrix::BlockType::const_iterator base_iterator;
	
					 /**
					  * Number of block where row lies in.
					  */
	unsigned int row_block;
	
					 /**
					  * Number of block column where
					  * column lies in.
					  */
	unsigned int col_block;
	
					 /**
					  * Index in whole row.
					  */
	unsigned int a_index;

                                         /**
                                          * Let the iterator class be a
                                          * friend.
                                          */
        template <typename>
        friend class ConstIterator;
    };
    
				     /**
				      * STL conforming iterator.
				      */
    template <class BlockMatrix>
    class ConstIterator
    {
      public:
                                         /**
                                          * Constructor.
                                          */ 
	ConstIterator(const BlockMatrix  *matrix,
                      const unsigned int  row,
                      const unsigned int  index);
	  
                                         /**
                                          * Prefix increment.
                                          */
	ConstIterator& operator++ ();

                                         /**
                                          * Postfix increment.
                                          */
	ConstIterator& operator++ (int);

                                         /**
                                          * Dereferencing operator.
                                          */
	const Accessor<BlockMatrix> & operator* () const;

                                         /**
                                          * Dereferencing operator.
                                          */
	const Accessor<BlockMatrix> * operator-> () const;

                                         /**
                                          * Comparison. True, if
                                          * both iterators point to
                                          * the same matrix
                                          * position.
                                          */
	bool operator == (const ConstIterator&) const;
                                         /**
                                          * Inverse of operator==().
                                          */
	bool operator != (const ConstIterator&) const;

                                         /**
                                          * Comparison
                                          * operator. Result is true
                                          * if either the first row
                                          * number is smaller or if
                                          * the row numbers are
                                          * equal and the first
                                          * index is smaller.
                                          */
	bool operator < (const ConstIterator&) const;

      private:
                                         /**
                                          * The accessor with which we work.
                                          */
        Accessor<BlockMatrix> accessor;
    };
  }
}


/**
 * Blocked matrix class. The behaviour of objects of this type is almost as
 * for the usual matrix objects, with most of the functions being implemented
 * in both classes. The main difference is that the matrix represented by this
 * object is composed of an array of matrices (e.g. of type
 * SparseMatrix<number>) and all accesses to the elements of this object are
 * relayed to accesses of the base matrices. The actual type of the individual
 * blocks of this matrix is the type of the template argument, and can, for
 * example be the usual SparseMatrix or PETScWrappers::SparseMatrix.
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
 * Objects of this type are frequently used when a system of differential
 * equations has solutions with variables that fall into different
 * classes. For example, solutions of the Stokes or Navier-Stokes equations
 * have @p dim velocity components and one pressure component. In this case,
 * it may make sense to consider the linear system of equations as a system of
 * 2x2 blocks, and one can construct preconditioners or solvers based on this
 * 2x2 block structure. This class can help you in these cases, as it allows
 * to view the matrix alternatively as one big matrix, or as a number of
 * individual blocks.
 *
 * 
 * @section Inheritance Inheriting from this class
 *
 * Since this class simply forwards its calls to the subobjects (if necessary
 * after adjusting indices denoting which subobject is meant), this class is
 * completely independent of the actual type of the subobject. The functions
 * that set up block matrices and destroy them, however, have to be
 * implemented in derived classes. These functions also have to fill the data
 * members provided by this base class, as they are only used passively in
 * this class.
 *
 * 
 * @ref Instantiations: some (<tt>@<float@> @<double@></tt>). Most of the
 * functions take a vector or block vector argument. These functions can, in
 * general, only successfully be compiled if the individual blocks of this
 * matrix implement the respective functions operating on the vector type in
 * question. For example, if you have a block sparse matrix over deal.II
 * SparseMatrix objects, then you will likely not be able to form the
 * matrix-vector multiplication with a block vector over
 * PETScWrappers::SparseMatrix objects. If you attempt anyway, you will likely
 * get a number of compiler errors.
 *
 * @author Wolfgang Bangerth, 2000, 2004
 */
template <typename MatrixType>
class BlockMatrixBase : public Subscriptor
{
  public:
                                     /**
                                      * Typedef the type of the underlying
                                      * matrix.
                                      */
    typedef MatrixType BlockType;

				     /**
				      * Type of matrix entries. In analogy to
				      * the STL container classes.
				      */
    typedef typename BlockType::value_type value_type;
    typedef value_type             *pointer;
    typedef const value_type       *const_pointer;
    typedef value_type             &reference;
    typedef const value_type       &const_reference;
    typedef size_t                  size_type;

    typedef
    internal::BlockMatrixIterators::ConstIterator<BlockMatrixBase>
    iterator;

    typedef
    internal::BlockMatrixIterators::ConstIterator<BlockMatrixBase>
    const_iterator;


                                     /**
                                      * Default constructor.
                                      */
    BlockMatrixBase ();
    
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
    template <class BlockMatrixType>
    BlockMatrixBase &
    copy_from (const BlockMatrixType &source);

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
				      * Access the block with the
				      * given coordinates.
				      */
    BlockType &
    block (const unsigned int row,
	   const unsigned int column);
    
    
				     /**
				      * Access the block with the
				      * given coordinates. Version for
				      * constant objects.
				      */
    const BlockType &
    block (const unsigned int row,
	   const unsigned int column) const;    
    
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
				      * Set the element <tt>(i,j)</tt>
				      * to <tt>value</tt>. Throws an
				      * error if the entry does not
				      * exist. Still, it is allowed to
				      * store zero values in
				      * non-existent fields.
				      */
    void set (const unsigned int i,
	      const unsigned int j,
	      const value_type value);
    
				     /**
				      * Add <tt>value</tt> to the element
				      * <tt>(i,j)</tt>.  Throws an error if
				      * the entry does not
				      * exist. Still, it is allowed to
				      * store zero values in
				      * non-existent fields.
				      */
    void add (const unsigned int i,
              const unsigned int j,
	      const value_type value);

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
    value_type operator () (const unsigned int i,
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
    value_type el (const unsigned int i,
                   const unsigned int j) const;

                                     /**
                                      * Call the compress() function on all
                                      * the subblocks of the matrix.
                                      */
    void compress ();

                                     /**
				      * Multiply the entire matrix by a
				      * fixed factor.
				      */
    BlockMatrixBase & operator *= (const value_type factor);
    
				     /**
				      * Divide the entire matrix by a
				      * fixed factor.
				      */
    BlockMatrixBase & operator /= (const value_type factor);
    
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
    template <class BlockMatrixType>
    void add_scaled (const value_type       factor,
		     const BlockMatrixType &matrix);
    
				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M*src$ with $M$
				      * being this matrix.
				      */
    template <class BlockVectorType>
    void vmult (BlockVectorType       &dst,
		const BlockVectorType &src) const;

				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block column.
				      */
    template <class BlockVectorType,
              class VectorType>
    void vmult (BlockVectorType          &dst,
		const VectorType &src) const;

				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block row.
				      */
    template <class BlockVectorType,
              class VectorType>
    void vmult (VectorType    &dst,
		const BlockVectorType &src) const;

				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block.
				      */
    template <class BlockVectorType,
              class VectorType>
    void vmult (VectorType       &dst,
		const VectorType &src) const;
    
				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M^T*src$ with $M$
				      * being this matrix. This
				      * function does the same as
				      * vmult() but takes the
				      * transposed matrix.
				      */
    template <class BlockVectorType>
    void Tvmult (BlockVectorType       &dst,
		 const BlockVectorType &src) const;
  
				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block row.
				      */
    template <class BlockVectorType,
              class VectorType>
    void Tvmult (BlockVectorType  &dst,
		 const VectorType &src) const;

				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block column.
				      */
    template <class BlockVectorType,
              class VectorType>
    void Tvmult (VectorType    &dst,
		 const BlockVectorType &src) const;

				     /**
				      * Matrix-vector
				      * multiplication. Just like the
				      * previous function, but only
				      * applicable if the matrix has
				      * only one block.
				      */
    template <class BlockVectorType,
              class VectorType>
    void Tvmult (VectorType       &dst,
		 const VectorType &src) const;
    
				     /**
				      * Adding Matrix-vector
				      * multiplication. Add $M*src$ on
				      * $dst$ with $M$ being this
				      * matrix.
				      */
    template <class BlockVectorType>
    void vmult_add (BlockVectorType       &dst,
		    const BlockVectorType &src) const;
    
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
    template <class BlockVectorType>
    void Tvmult_add (BlockVectorType       &dst,
		     const BlockVectorType &src) const;
  
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
    template <class BlockVectorType>
    value_type
    matrix_norm_square (const BlockVectorType &v) const;

				     /**
				      * Compute the matrix scalar
				      * product $\left(u,Mv\right)$.
				      */    
    template <class BlockVectorType>
    value_type
    matrix_scalar_product (const BlockVectorType &u,
			   const BlockVectorType &v) const;
    
				     /**
				      * Compute the residual
				      * <i>r=b-Ax</i>. Write the
				      * residual into <tt>dst</tt>.
				      */
    template <class BlockVectorType>
    value_type residual (BlockVectorType       &dst,
			 const BlockVectorType &x,
			 const BlockVectorType &b) const;

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
    const_iterator begin (const unsigned int r) const;

				     /**
				      * Final iterator of row <tt>r</tt>.
				      */
    const_iterator end (const unsigned int r) const;


				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotBlockSquare);

				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotInitialized);
                                     /**
                                      * Exception
                                      */
    DeclException4 (ExcIncompatibleRowNumbers,
                    int, int, int, int,
                    << "The blocks [" << arg1 << ',' << arg2 << "] and ["
                    << arg3 << ',' << arg4 << "] have differing row numbers.");
                                     /**
                                      * Exception
                                      */
    DeclException4 (ExcIncompatibleColNumbers,
                    int, int, int, int,
                    << "The blocks [" << arg1 << ',' << arg2 << "] and ["
                    << arg3 << ',' << arg4 << "] have differing column numbers.");

  protected:
                                     /**
                                      * Index arrays for rows and columns.
                                      */
    BlockIndices row_block_indices;
    BlockIndices column_block_indices;
    
				     /**
				      * Array of sub-matrices.
				      */
    Table<2,SmartPointer<BlockType> > sub_objects;

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
                                      *
                                      * Derived classes should call this
                                      * function whenever the size of the
                                      * sub-objects has changed and the @p
                                      * X_block_indices arrays need to be
                                      * updated.
                                      *
                                      * Note that this function is not public
                                      * since not all derived classes need to
                                      * export its interface. For example, for
                                      * the usual deal.II SparseMatrix class,
                                      * the sizes are implicitly determined
                                      * whenever reinit() is called, and
                                      * individual blocks cannot be
                                      * resized. For that class, this function
                                      * therefore does not have to be
                                      * public. On the other hand, for the
                                      * PETSc classes, there is no associated
                                      * sparsity pattern object that
                                      * determines the block sizes, and for
                                      * these the function needs to be
                                      * publicly available. These classes
                                      * therefore export this function.
                                      */
    void collect_sizes ();
    
				     /**
				      * Make the iterator class a
				      * friend. We have to work around
				      * a compiler bug here again.
				      */
#ifndef DEAL_II_NAMESP_TEMPL_FRIEND_BUG
    template <typename>
    friend class internal::BlockMatrixIterators::Accessor;

    template <typename>
    friend class internal::BlockMatrixIterators::ConstIterator;
#else
    typedef internal::BlockMatrixIterators::Accessor<BlockMatrixBase> Accessor;
    friend class Accessor;
    
    friend class const_iterator;
#endif    
};


/*@}*/
/* ------------------------- Template functions ---------------------- */


namespace internal
{
  namespace BlockMatrixIterators
  {
    template <class BlockMatrix>
    inline
    Accessor<BlockMatrix>::
    Accessor (const BlockMatrix  *matrix,
              const unsigned int  r,
              const unsigned int  i)
                    :
                    matrix(matrix),
                    base_iterator(matrix->block(0,0).begin()),
                    row_block(0),
                    col_block(0),
                    a_index(0)
    {
      Assert (i==0, ExcNotImplemented());

      if (r < matrix->m())
        {
          std::pair<unsigned int,unsigned int> indices
            = matrix->row_block_indices.global_to_local(r);
          row_block = indices.first;
          base_iterator = matrix->block(indices.first, 0).begin(indices.second);
        }
      else
        {
          row_block = matrix->n_block_rows();
          base_iterator = matrix->block(0, 0).begin();
        }
    }


    template <class BlockMatrix>
    inline
    unsigned int
    Accessor<BlockMatrix>::row() const
    {
      if (row_block < matrix->n_block_rows())
        return (matrix->row_block_indices.local_to_global(row_block, 0) +
                base_iterator->row());
      else
        return 0;
    }


    template <class BlockMatrix>
    inline
    unsigned int
    Accessor<BlockMatrix>::index() const
    {
      return a_index;
    }


    template <class BlockMatrix>
    inline
    unsigned int
    Accessor<BlockMatrix>::column() const
    {
                                       // note: the block column should always
                                       // be valid, in contrast to the block
                                       // row, which may be past the end of
                                       // the matrix. thus, unlike in the
                                       // row() function, we do not have to
                                       // if-else here
      return (matrix->column_block_indices.local_to_global(col_block,0) +
              base_iterator->column());
    }


    template <class BlockMatrix>
    inline
    unsigned int
    Accessor<BlockMatrix>::block_row() const
    {
      return row_block;
    }


    template <class BlockMatrix>
    inline
    unsigned int
    Accessor<BlockMatrix>::block_column() const
    {
      return col_block;
    }


    template <class BlockMatrix>
    inline
    typename Accessor<BlockMatrix>::value_type
    Accessor<BlockMatrix>::value () const
    {
      return base_iterator->value();
    }


//----------------------------------------------------------------------//


    template <class BlockMatrix>
    inline
    ConstIterator<BlockMatrix>::
    ConstIterator (const BlockMatrix *m,
                   const unsigned int r,
                   const unsigned int i)
                    :
                    accessor (m, r, i)
    {}



    template <class BlockMatrix>
    inline
    ConstIterator<BlockMatrix> &
    ConstIterator<BlockMatrix>::operator++ ()
    {
      Assert (accessor.row_block < accessor.matrix->n_block_rows(),
              ExcIteratorPastEnd());

                                       // Remember current row inside block
      unsigned int local_row = accessor.base_iterator->row();

                                       // Advance inside block
      ++accessor.base_iterator;
      ++accessor.a_index;
  
                                       // If end of row inside block,
                                       // advance to next block
      if (accessor.base_iterator ==
          accessor.matrix->block(accessor.row_block,
                                 accessor.col_block).end(local_row))
        {
          if (accessor.col_block < accessor.matrix->n_block_cols()-1)
            {
                                               // Advance to next block in
                                               // row
//TODO: what if this row in that block is empty?          
              ++accessor.col_block;
            }
          else
            {
                                               // Advance to first block
                                               // in next row
              accessor.col_block = 0;
              accessor.a_index = 0;
              ++local_row;
              if (local_row >=
                  accessor.matrix->block(accessor.row_block,0).m())
                {
                                                   // If final row in
                                                   // block, go to next
                                                   // block row
                  local_row = 0;
                  ++accessor.row_block;
                }
            }
                                           // Finally, set base_iterator
                                           // to start of row determined
                                           // above
          if (accessor.row_block < accessor.matrix->n_block_rows())
            accessor.base_iterator =
              accessor.matrix->block(accessor.row_block, accessor.col_block)
              .begin(local_row);
          else
                                             // Set base_iterator to a
                                             // defined state for
                                             // comparison. This is the
                                             // end() state.
//TODO: this is a particularly bad choice, since it is actually a valid iterator. it should rather be something like the end iterator of the last block!
            accessor.base_iterator = accessor.matrix->block(0, 0).begin();
        }
      return *this;
    }



//  template <class BlockMatrix>
//  inline
//  const_iterator&
//  ConstIterator::operator++ (int)
//  {
//    Assert (false, ExcNotImplemented());
//  }




    template <class BlockMatrix>
    inline
    const Accessor<BlockMatrix> &
    ConstIterator<BlockMatrix>::operator* () const
    {
      return accessor;
    }


    template <class BlockMatrix>
    inline
    const Accessor<BlockMatrix> *
    ConstIterator<BlockMatrix>::operator-> () const
    {
      return &accessor;
    }



    template <class BlockMatrix>
    inline
    bool
    ConstIterator<BlockMatrix>::
    operator == (const ConstIterator& i) const
    {
      if (accessor.matrix != i.accessor.matrix)
        return false;
  
      if (accessor.row_block == i.accessor.row_block
          && accessor.col_block == i.accessor.col_block
          && accessor.base_iterator == i.accessor.base_iterator)
        return true;
      return false;
    }



    template <class BlockMatrix>
    inline
    bool
    ConstIterator<BlockMatrix>::
    operator != (const ConstIterator& i) const
    {
      return !(*this == i);
    }



    template <class BlockMatrix>
    inline
    bool
    ConstIterator<BlockMatrix>::
    operator < (const ConstIterator& i) const
    {
      if (accessor.row_block < i.accessor.row_block)
        return true;
      if (accessor.row_block == i.accessor.row_block)
        {
          if (accessor.base_iterator->row() < i.accessor.base_iterator->row())
            return true;
          if (accessor.base_iterator->row() == i.accessor.base_iterator->row())
            if (accessor.a_index < i.accessor.a_index)
              return true;
        }
      return false;
    }

  }
}

//----------------------------------------------------------------------//


template <typename MatrixType>
inline
BlockMatrixBase<MatrixType>::BlockMatrixBase ()
{}



template <class MatrixType>
template <class BlockMatrixType>
inline
BlockMatrixBase<MatrixType> &
BlockMatrixBase<MatrixType>::
copy_from (const BlockMatrixType &source)
{
  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c).copy_from (source.block(r,c));
  
  return *this;
}



template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::clear () 
{
  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c).clear ();
  sub_objects.reinit (0,0);
  row_block_indices = column_block_indices = BlockIndices();

                                   // reset block indices to empty
  row_block_indices = column_block_indices = BlockIndices ();
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::BlockType &
BlockMatrixBase<MatrixType>::block (const unsigned int row,
                                    const unsigned int column)
{
  Assert (row<n_block_rows(),
          ExcIndexRange (row, 0, n_block_rows()));
  Assert (column<n_block_cols(),
          ExcIndexRange (column, 0, n_block_cols()));
  
  return *sub_objects[row][column];
}



template <class MatrixType>
inline
const typename BlockMatrixBase<MatrixType>::BlockType &
BlockMatrixBase<MatrixType>::block (const unsigned int row,
                                   const unsigned int column) const
{
  Assert (row<n_block_rows(),
          ExcIndexRange (row, 0, n_block_rows()));
  Assert (column<n_block_cols(),
          ExcIndexRange (column, 0, n_block_cols()));
  
  return *sub_objects[row][column];
}


template <class MatrixType>
inline
unsigned int
BlockMatrixBase<MatrixType>::m () const
{
  return row_block_indices.total_size();
}



template <class MatrixType>
inline
unsigned int
BlockMatrixBase<MatrixType>::n () const
{
  return column_block_indices.total_size();
}



template <class MatrixType>
inline
unsigned int
BlockMatrixBase<MatrixType>::n_block_cols () const
{
  return column_block_indices.size();
}



template <class MatrixType>
inline
unsigned int
BlockMatrixBase<MatrixType>::n_block_rows () const
{
  return row_block_indices.size();
}



template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::set (const unsigned int i,
                                  const unsigned int j,
                                  const value_type value)
{
  const std::pair<unsigned int,unsigned int>
    row_index = row_block_indices.global_to_local (i),
    col_index = column_block_indices.global_to_local (j);
  block(row_index.first,col_index.first).set (row_index.second,
					      col_index.second,
					      value);
}



template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::add (const unsigned int i,
                                  const unsigned int j,
                                  const value_type value)
{
                                   // save some cycles for zero additions, but
                                   // only if it is safe for the matrix we are
                                   // working with
  if ((MatrixType::Traits::zero_addition_can_be_elided == true)
      &&
      (value == 0))
    return;
  
  const std::pair<unsigned int,unsigned int>
    row_index = row_block_indices.global_to_local (i),
    col_index = column_block_indices.global_to_local (j);
  block(row_index.first,col_index.first).add (row_index.second,
					      col_index.second,
					      value);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::operator () (const unsigned int i,
                                          const unsigned int j) const
{
  const std::pair<unsigned int,unsigned int>
    row_index = row_block_indices.global_to_local (i),
    col_index = column_block_indices.global_to_local (j);
  return block(row_index.first,col_index.first) (row_index.second,
						 col_index.second);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::el (const unsigned int i,
                                 const unsigned int j) const
{
  const std::pair<unsigned int,unsigned int>
    row_index = row_block_indices.global_to_local (i),
    col_index = column_block_indices.global_to_local (j);
  return block(row_index.first,col_index.first).el (row_index.second,
						    col_index.second);
}



template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::compress ()
{
  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c).compress ();
}



template <class MatrixType>
inline
BlockMatrixBase<MatrixType> &
BlockMatrixBase<MatrixType>::operator *= (const value_type factor)
{
  Assert (n_block_cols() != 0, ExcMatrixNotInitialized());
  Assert (n_block_rows() != 0, ExcMatrixNotInitialized());

  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c) *= factor;

  return *this;
}



template <class MatrixType>
inline
BlockMatrixBase<MatrixType> &
BlockMatrixBase<MatrixType>::operator /= (const value_type factor)
{
  Assert (n_block_cols() != 0, ExcMatrixNotInitialized());
  Assert (n_block_rows() != 0, ExcMatrixNotInitialized());
  Assert (factor !=0, ExcDivideByZero());

  const value_type factor_inv = 1. / factor;

  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c) *= factor_inv;

  return *this;
}



template <class MatrixType>
template <class BlockMatrixType>
void
BlockMatrixBase<MatrixType>::
add_scaled (const value_type factor,
	    const BlockMatrixType &matrix)
{
  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c).add_scaled (factor, matrix.block(r,c));
}




template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::vmult (BlockVectorType       &dst,
                                    const BlockVectorType &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  for (unsigned int row=0; row<n_block_rows(); ++row)
    {
      block(row,0).vmult (dst.block(row),
			  src.block(0));
      for (unsigned int col=1; col<n_block_cols(); ++col)
	block(row,col).vmult_add (dst.block(row),
				  src.block(col));
    };
}



template <class MatrixType>
template <class BlockVectorType,
          class VectorType>
void
BlockMatrixBase<MatrixType>::vmult (VectorType    &dst,
                                    const BlockVectorType &src) const
{
  Assert (n_block_rows() == 1,
	  ExcDimensionMismatch(1, n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  block(0,0).vmult (dst, src.block(0));
  for (unsigned int col=1; col<n_block_cols(); ++col)
    block(0,col).vmult_add (dst, src.block(col));
}



template <class MatrixType>
template <class BlockVectorType,
          class VectorType>
void
BlockMatrixBase<MatrixType>::vmult (BlockVectorType  &dst,
                                    const VectorType &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (1 == n_block_cols(),
	  ExcDimensionMismatch(1, n_block_cols()));

  for (unsigned int row=0; row<n_block_rows(); ++row)
    block(row,0).vmult (dst.block(row),
			src);
}



template <class MatrixType>
template <class BlockVectorType,
          class VectorType>
void
BlockMatrixBase<MatrixType>::vmult (VectorType       &dst,
                                    const VectorType &src) const
{
  Assert (1 == n_block_rows(),
	  ExcDimensionMismatch(1, n_block_rows()));
  Assert (1 == n_block_cols(),
	  ExcDimensionMismatch(1, n_block_cols()));

  block(0,0).vmult (dst, src);
}



template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::vmult_add (BlockVectorType       &dst,
                                        const BlockVectorType &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  for (unsigned int row=0; row<n_block_rows(); ++row)
    {
      block(row,0).vmult_add (dst.block(row),
			      src.block(0));
      for (unsigned int col=1; col<n_block_cols(); ++col)
	block(row,col).vmult_add (dst.block(row),
				  src.block(col));
    };
}




template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::Tvmult (BlockVectorType       &dst,
                                     const BlockVectorType &src) const
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



template <class MatrixType>
template <class BlockVectorType,
          class VectorType>
void
BlockMatrixBase<MatrixType>::Tvmult (BlockVectorType  &dst,
                                     const VectorType &src) const
{
  Assert (dst.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert (1 == n_block_rows(),
	  ExcDimensionMismatch(1, n_block_rows()));

  dst = 0.;
  
  for (unsigned int col=0; col<n_block_cols(); ++col)
    block(0,col).Tvmult_add (dst.block(col), src);
}



template <class MatrixType>
template <class BlockVectorType,
          class VectorType>
void
BlockMatrixBase<MatrixType>::Tvmult (VectorType    &dst,
                                     const BlockVectorType &src) const
{
  Assert (1 == n_block_cols(),
	  ExcDimensionMismatch(1, n_block_cols()));
  Assert (src.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  block(0,0).Tvmult (dst, src.block(0));
  
  for (unsigned int row=1; row<n_block_rows(); ++row)
    block(row,0).Tvmult_add (dst, src.block(row));
}



template <class MatrixType>
template <class BlockVectorType,
          class VectorType>
void
BlockMatrixBase<MatrixType>::Tvmult (VectorType       &dst,
                                     const VectorType &src) const
{
  Assert (1 == n_block_cols(),
	  ExcDimensionMismatch(1, n_block_cols()));
  Assert (1 == n_block_rows(),
	  ExcDimensionMismatch(1, n_block_rows()));

  block(0,0).Tvmult (dst, src);
}



template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::Tvmult_add (BlockVectorType& dst,
                                         const BlockVectorType& src) const
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



template <class MatrixType>
template <class BlockVectorType>
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::matrix_norm_square (const BlockVectorType &v) const
{
  Assert (n_block_rows() == n_block_cols(),
          ExcMatrixNotBlockSquare());
  Assert (v.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(v.n_blocks(), n_block_rows()));

  value_type norm_sqr = 0;
  for (unsigned int row=0; row<n_block_rows(); ++row)
    for (unsigned int col=0; col<n_block_cols(); ++col)
      if (row==col)
	norm_sqr += block(row,col).matrix_norm_square (v.block(row));
      else
	norm_sqr += block(row,col).matrix_scalar_product (v.block(row),
							  v.block(col));
  return norm_sqr;
}



template <class MatrixType>
template <class BlockVectorType>
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::
matrix_scalar_product (const BlockVectorType    &u,
		       const BlockVectorType &v) const
{
  Assert (u.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(u.n_blocks(), n_block_rows()));
  Assert (v.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(v.n_blocks(), n_block_cols()));

  value_type result = 0;
  for (unsigned int row=0; row<n_block_rows(); ++row)
    for (unsigned int col=0; col<n_block_cols(); ++col)
      result += block(row,col).matrix_scalar_product (u.block(row),
						      v.block(col));
  return result;
}



template <class MatrixType>
template <class BlockVectorType>
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::
residual (BlockVectorType          &dst,
	  const BlockVectorType &x,
	  const BlockVectorType    &b) const
{
  Assert (dst.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (b.n_blocks() == n_block_rows(),
	  ExcDimensionMismatch(b.n_blocks(), n_block_rows()));
  Assert (x.n_blocks() == n_block_cols(),
	  ExcDimensionMismatch(x.n_blocks(), n_block_cols()));
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
  for (unsigned int row=0; row<n_block_rows(); ++row)
    {
      block(row,0).residual (dst.block(row),
			     x.block(0),
			     b.block(row));
      
      for (unsigned int i=0; i<dst.block(row).size(); ++i)
	dst.block(row)(i) = -dst.block(row)(i);
      
      for (unsigned int col=1; col<n_block_cols(); ++col)
	block(row,col).vmult_add (dst.block(row),
				  x.block(col));

      for (unsigned int i=0; i<dst.block(row).size(); ++i)
	dst.block(row)(i) = -dst.block(row)(i);
    };

  value_type res = 0;
  for (unsigned int row=0; row<n_block_rows(); ++row)
    res += dst.block(row).norm_sqr ();
  return std::sqrt(res);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::begin () const
{
  return const_iterator(this, 0, 0);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::end () const
{
  return const_iterator(this, m(), 0);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::begin (const unsigned int r) const
{
  Assert (r<m(), ExcIndexRange(r,0,m()));
  return const_iterator(this, r, 0);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::end (const unsigned int r) const
{
  Assert (r<m(), ExcIndexRange(r,0,m()));
  return const_iterator(this, r+1, 0);
}



template <class MatrixType>
void
BlockMatrixBase<MatrixType>::collect_sizes ()
{
  std::vector<unsigned int> row_sizes (this->n_block_rows());
  std::vector<unsigned int> col_sizes (this->n_block_cols());

                                   // first find out the row sizes
                                   // from the first block column
  for (unsigned int r=0; r<this->n_block_rows(); ++r)
    row_sizes[r] = sub_objects[r][0]->m();
                                   // then check that the following
                                   // block columns have the same
                                   // sizes
  for (unsigned int c=1; c<this->n_block_cols(); ++c)
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      Assert (row_sizes[r] == sub_objects[r][c]->m(),
              ExcIncompatibleRowNumbers (r,0,r,c));

                                   // finally initialize the row
                                   // indices with this array
  this->row_block_indices.reinit (row_sizes);
  
  
                                   // then do the same with the columns
  for (unsigned int c=0; c<this->n_block_cols(); ++c)
    col_sizes[c] = sub_objects[0][c]->n();
  for (unsigned int r=1; r<this->n_block_rows(); ++r)
    for (unsigned int c=0; c<this->n_block_cols(); ++c)
      Assert (col_sizes[c] == sub_objects[r][c]->n(),
              ExcIncompatibleRowNumbers (0,c,r,c));

                                   // finally initialize the row
                                   // indices with this array
  this->column_block_indices.reinit (col_sizes);
}



#endif    // __deal2__block_matrix_base_h
