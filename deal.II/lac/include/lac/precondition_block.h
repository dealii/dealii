//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------
#ifndef __deal2__precondition_block_h
#define __deal2__precondition_block_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>

#include <vector>

template <typename number> class FullMatrix;
template <typename number> class Vector;

template<class MATRIX, typename inverse_type>
class PreconditionBlockJacobi;

/*! @addtogroup Preconditioners
 *@{
 */


/**
 * Base class for @p PreconditionBlockJacobi,
 * @p PreconditionBlockSOR, ...  This class assumes the
 * <tt>SparseMatrix<number></tt> consisting of invertible blocks of
 * @p blocksize on the diagonal and provides the inversion of the
 * diagonal blocks of the matrix. NOT only block diagonal matrices are
 * allowed but all matrices of arbitrary structure with the minimal
 * property of having invertible blocks on the diagonal!
 *
 * This block matrix structure is given e.g. for the DG method for the
 * transport equation. For a downstream numbering the matrices even
 * have got a block lower left matrix structure, i.e. the matrices are
 * empty above the diagonal blocks.
 *
 * For all matrices that are empty above and below the diagonal blocks
 * (i.e. for all block diagonal matrices) the @p BlockJacobi
 * preconditioner is a direct solver. For all matrices that are empty
 * only above the diagonal blocks (e.g. the matrices one gets by the
 * DG method with downstream numbering) @p BlockSOR is a direct
 * solver.
 * 
 * This first implementation of the @p PreconditionBlock assumes the
 * matrix has blocks each of the same block size. Varying block sizes
 * within the matrix must still be implemented if needed.
 *
 * The first template parameter denotes the type of number
 * representation in the sparse matrix, the second denotes the type of
 * number representation in which the inverted diagonal block matrices
 * are stored within this class by <tt>invert_diagblocks()</tt>. If you
 * don't want to use the block inversion as an exact solver, but
 * rather as a preconditioner, you may probably want to store the
 * inverted blocks with less accuracy than the original matrix; for
 * example, <tt>number==double, inverse_type=float</tt> might be a viable
 * choice.
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
 * @author Ralf Hartmann, Guido Kanschat, 1999, 2000
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class PreconditionBlock : public virtual Subscriptor
{
  private:
				     /**
				      * Define number type of matrix.
				      */
    typedef typename MATRIX::value_type number;

				     /**
				      * Value type for inverse matrices.
				      */
    typedef inverse_type value_type;
    
  public:
				     /**
				      * Parameters for block preconditioners.
				      */
    class AdditionalData
    {
      public:
					 /**
					  * Constructor. Block size
					  * must be given since there
					  * is no reasonable default
					  * parameter.
					  */
	AdditionalData (const unsigned int block_size,
			const double relaxation = 1.,
			const bool invert_diagonal = true,
			const bool same_diagonal = false);

					 /**
					  * Relaxation parameter.
					  */
	double relaxation;
	
					 /**
					  * Block size.
					  */
	unsigned int block_size;

					 /**
					  * Invert diagonal during initialization.
					  */
	bool invert_diagonal;

					 /**
					  * Assume all diagonal blocks
					  * are equal to save memory.
					  */
	bool same_diagonal;
    };
    
    
				     /**
				      * Constructor.
				      */
    PreconditionBlock();
    
				     /**
				      * Destructor.
				      */
    ~PreconditionBlock();

				     /**
				      * Initialize matrix and block
				      * size.  We store the matrix and
				      * the block size in the
				      * preconditioner object. In a
				      * second step, the inverses of
				      * the diagonal blocks may be
				      * computed.
				      *
				      * Additionally, a relaxation
				      * parameter for derived classes
				      * may be provided.
				      */
    void initialize (const MATRIX& A,
		     const AdditionalData parameters);


				     /**
				      * Initialize matrix and block
				      * size for permuted
				      * preconditioning. Additionally
				      * to the parameters of the other
				      * initalize() function, we hand
				      * over two index vectors with
				      * the permutation and its
				      * inverse. For the meaning of
				      * these vectors see
				      * PreconditionBlockSOR.
				      *
				      * In a second step, the inverses
				      * of the diagonal blocks may be
				      * computed. Make sure you use
				      * invert_permuted_diagblocks()
				      * to yield consistent data.
				      *
				      * Additionally, a relaxation
				      * parameter for derived classes
				      * may be provided.
				      */
    void initialize (const MATRIX& A,
		     const std::vector<unsigned int>& permutation,
		     const std::vector<unsigned int>& inverse_permutation,
		     const AdditionalData parameters);

				     /**
				      * Replacement of
				      * invert_diagblocks() for
				      * permuted preconditioning.
				      */ 
    void invert_permuted_diagblocks(
      const std::vector<unsigned int>& permutation,
      const std::vector<unsigned int>& inverse_permutation);

				     /**
				      * Deletes the inverse diagonal
				      * block matrices if existent,
				      * sets the blocksize to 0, hence
				      * leaves the class in the state
				      * that it had directly after
				      * calling the constructor.
				      */
    void clear();

				     /**
				      * Checks whether the object is empty.
				      */
    bool empty () const;

				     /**
				      * Read-only access to entries.
				      * This function is only possible
				      * if the inverse diagonal blocks
				      * are stored.
				      */
    value_type el(unsigned int i,
		  unsigned int j) const;
    
				     /**
				      * Use only the inverse of the
				      * first diagonal block to save
				      * memory and computation time.
				      *
				      * Possible applications:
				      * computing on a cartesian grid,
				      * all diagonal blocks are the
				      * same or all diagonal blocks
				      * are at least similar and
				      * inversion of one of them still
				      * yields a preconditioner.
				      */
    void set_same_diagonal ();

				     /**
				      * Does the matrix use only one
				      * diagonal block?
				      */
    bool same_diagonal () const;
    
    				     /**
				      * Stores the inverse of the
				      * diagonal blocks in
				      * @p inverse. This costs some
				      * additional memory - for DG
				      * methods about 1/3 (for double
				      * inverses) or 1/6 (for float
				      * inverses) of that used for the
				      * matrix - but it makes the
				      * preconditioning much faster.
				      *
				      * It is not allowed to call this
				      * function twice (will produce
				      * an error) before a call of
				      * <tt>clear(...)</tt>  because at the
				      * second time there already
				      * exist the inverse matrices.
				      *
				      * After this function is called,
				      * the lock on the matrix given
				      * through the @p use_matrix
				      * function is released, i.e. you
				      * may overwrite of delete it.
				      * You may want to do this in
				      * case you use this matrix to
				      * precondition another matrix.
				      */
    void invert_diagblocks();

				     /**
				      * Return the size of the blocks.
				      */
    unsigned int block_size () const;

				     /**
				      * The number of blocks of the
				      * matrix.
				      */
    unsigned int n_blocks() const;
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * Determine, whether inverses
				      * have been computed.
				      */
    bool inverses_ready () const;
    
				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongBlockSize,
		    int, int,
		    << "The blocksize " << arg1
		    << " and the size of the matrix " << arg2
		    << " do not match.");

				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongNumberOfInverses,
		    int, int,
		    << "There are " << arg1
		    << " inverse matrices but " << arg2
		    << " cells.");

				     /**
				      * Exception
				      */
    DeclException0 (ExcInverseMatricesAlreadyExist);

				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotSquare);

				     /**
				      * Exception
				      */
    DeclException0 (ExcDiagonalsNotStored);
   
				     /**
				      * Access to the inverse diagonal
				      * blocks.
				      */
    const FullMatrix<inverse_type>& inverse (unsigned int i) const;
    
				     /**
				      * Access to the diagonal
				      * blocks.
				      */
    const FullMatrix<inverse_type>& diagonal (unsigned int i) const;
    
  protected:
				     /**
				      * Size of the blocks. Each
				      * diagonal block is assumed to
				      * be of the same size.
				      */
    unsigned int blocksize;

				     /**
				      * Pointer to the matrix. Make
				      * sure that the matrix exists as
				      * long as this class needs it,
				      * i.e. until calling
				      * @p invert_diagblocks, or (if
				      * the inverse matrices should
				      * not be stored) until the last
				      * call of the preconditoining
				      * @p vmult function of the
				      * derived classes.
				      */
    SmartPointer<const MATRIX> A;
				     /**
				      * Relaxation parameter to be
				      * used by derived classes.
				      */
    double relaxation;


				     /**
				      * Flag for storing the diagonal
				      * blocks of the matrix.
				      */
    bool store_diagonals;

				     /**
				      * Number of blocks.
				      */
    unsigned int nblocks;
  private:
    
				     /**
				      * Storage of the inverse
				      * matrices of the diagonal
				      * blocks matrices as
				      * <tt>FullMatrix<inverse_type></tt>
				      * matrices. Using
				      * <tt>inverse_type=float</tt> saves
				      * memory in comparison with
				      * <tt>inverse_type=double</tt>.
				      */
    std::vector<FullMatrix<inverse_type> > var_inverse;

				     /**
				      * Storage of the original diagonal blocks.
				      * These are only filled if @p store_diagonals
				      * is @p true.
				      *
				      * Used by the blocked SSOR method.
				      */
    std::vector<FullMatrix<inverse_type> > var_diagonal;
				      
				     /**
				      * Flag for diagonal compression.
				      * @ref set_same_diagonal()
				      */
    bool var_same_diagonal;
};



/**
 * Block Jacobi preconditioning.
 * See PreconditionBlock for requirements on the matrix.
 * @author Ralf Hartmann, Guido Kanschat, 1999, 2000, 2003
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class PreconditionBlockJacobi : public virtual Subscriptor,
				private PreconditionBlock<MATRIX, inverse_type>
{
  private:
				     /**
				      * Define number type of matrix.
				      */
    typedef typename MATRIX::value_type number;
    
  public:
				     /**
				      * STL conforming iterator.
				      */
    class const_iterator
    {
      private:
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
            Accessor (const PreconditionBlockJacobi<MATRIX, inverse_type> *matrix,
                      const unsigned int row);

                                             /**
                                              * Row number of the element
                                              * represented by this
                                              * object.
                                              */
            unsigned int row() const;

                                             /**
                                              * Column number of the
                                              * element represented by
                                              * this object.
                                              */
            unsigned int column() const;

                                             /**
                                              * Value of this matrix entry.
                                              */
            inverse_type value() const;
	
          protected:
                                             /**
                                              * The matrix accessed.
                                              */
            const PreconditionBlockJacobi<MATRIX, inverse_type>* matrix;

					     /**
					      * Save block size here
					      * for further reference.
					      */
	    unsigned int bs;
	    
					     /**
					      * Current block number.
					      */
	    unsigned int a_block;

					     /**
					      * Iterator inside block.
					      */
	    typename FullMatrix<inverse_type>::const_iterator b_iterator;

					     /**
					      * End of current block.
					      */
	    typename FullMatrix<inverse_type>::const_iterator b_end;
	    
                                             /**
                                              * Make enclosing class a
                                              * friend.
                                              */
            friend class const_iterator;
        };
        
      public:
                                         /**
                                          * Constructor.
                                          */ 
	const_iterator(const PreconditionBlockJacobi<MATRIX, inverse_type>* matrix,
		       const unsigned int row);
	  
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
                                          * Inverse of <tt>==</tt>.
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

      private:
                                         /**
                                          * Store an object of the
                                          * accessor class.
                                          */
        Accessor accessor;
    };

				     /**
				      * Make initialization function
				      * publicly available.
				      */
    PreconditionBlock<MATRIX, inverse_type>::initialize;
    
				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<MATRIX, inverse_type>::clear;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<MATRIX, inverse_type>::empty;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<MATRIX, inverse_type>::el;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<MATRIX, inverse_type>::set_same_diagonal;

				     /**
				      * Make function public.
				      */
    PreconditionBlock<MATRIX, inverse_type>::invert_diagblocks;
    
				     /**
				      * Make function public.
				      */
    PreconditionBlock<MATRIX, inverse_type>::block_size;
    
				     /**
				      * Make function public.
				      */
    PreconditionBlock<MATRIX, inverse_type>::n_blocks;
				     /**
				      * Make function accessible to iterator.
				      */
    PreconditionBlock<MATRIX, inverse_type>::inverse;
    
				     /**
				      * Execute block Jacobi
				      * preconditioning.
				      *
				      * This function will
				      * automatically use the inverse
				      * matrices if they exist, if not
				      * then BlockJacobi will need
				      * much time inverting the
				      * diagonal block matrices in
				      * each preconditioning step.
				      */
    template <typename number2>
    void vmult (Vector<number2>&, const Vector<number2>&) const;

				     /**
				      * Same as @p vmult, since Jacobi is symmetric.
				      */
    template <typename number2>
    void Tvmult (Vector<number2>&, const Vector<number2>&) const;
				     /**
				      * Execute block Jacobi
				      * preconditioning, adding to @p dst.
				      *
				      * This function will
				      * automatically use the inverse
				      * matrices if they exist, if not
				      * then BlockJacobi will need
				      * much time inverting the
				      * diagonal block matrices in
				      * each preconditioning step.
				      */
    template <typename number2>
    void vmult_add (Vector<number2>&, const Vector<number2>&) const;

				     /**
				      * Same as @p vmult_add, since Jacobi is symmetric.
				      */
    template <typename number2>
    void Tvmult_add (Vector<number2>&, const Vector<number2>&) const;

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
				      * first entry of row @p r.
				      */
    const_iterator begin (const unsigned int r) const;

				     /**
				      * Final iterator of row @p r.
				      */
    const_iterator end (const unsigned int r) const;
    

  private:
				   /**
				    * Actual implementation of the
				    * preconditioner.
				    *
				    * Depending on @p adding, the
				    * result of preconditioning is
				    * added to the destination vector.
				    */
    template <typename number2>
    void do_vmult (Vector<number2>&,
		   const Vector<number2>&,
		   bool adding) const;

    friend class Accessor;
    friend class const_iterator;
};



/**
 * Block SOR preconditioning.
 *
 * The functions @p vmult and @p Tvmult execute a (transposed)
 * block-SOR step, based on the blocks in PreconditionBlock. The
 * elements outside the diagonal blocks may be distributed
 * arbitrarily.
 *
 * See PreconditionBlock for requirements on the matrix.
 * @author Ralf Hartmann, Guido Kanschat, 1999, 2000, 2001, 2002, 2003
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class PreconditionBlockSOR : public virtual Subscriptor,
			     protected PreconditionBlock<MATRIX, inverse_type>
{
  private:
				     /**
				      * Define number type of matrix.
				      */
    typedef typename MATRIX::value_type number;
    
  public:
				     /**
				      * Make initialization function
				      * publicly available.
				      */
    PreconditionBlock<MATRIX, inverse_type>::initialize;
    
				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<MATRIX, inverse_type>::clear;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<MATRIX, inverse_type>::empty;
    
				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<MATRIX, inverse_type>::el;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<MATRIX, inverse_type>::set_same_diagonal;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<MATRIX, inverse_type>::invert_diagblocks;

				     /**
				      * Execute block SOR
				      * preconditioning.
				      *
				      * This function will
				      * automatically use the inverse
				      * matrices if they exist, if not
				      * then BlockSOR will waste much
				      * time inverting the diagonal
				      * block matrices in each
				      * preconditioning step.
				      *
				      * For matrices which are empty
				      * above the diagonal blocks
				      * BlockSOR is a direct solver.
				      */
    template <typename number2>
    void vmult (Vector<number2>&, const Vector<number2>&) const;

				     /**
				      * Execute block SOR
				      * preconditioning.
				      *
				      * Warning: this function
				      * performs normal @p vmult
				      * without adding. The reason for
				      * its existence is that
				      * BlockMatrixArray
				      * requires the adding version by
				      * default. On the other hand,
				      * adding requires an additional
				      * auxiliary vector, which is not
				      * desirable.
				      *
				      * @see vmult
				      */
    template <typename number2>
    void vmult_add (Vector<number2>&, const Vector<number2>&) const;

				     /**
				      * Backward application of vmult().
				      *
				      * In the current implementation,
				      * this is not the transpose of
				      * vmult(). It is a
				      * transposed Gauss-Seidel
				      * algorithm applied to the whole
				      * matrix, but the diagonal
				      * blocks being inverted are not
				      * transposed. Therefore, it is
				      * the transposed, if the
				      * diagonal blocks are symmetric.
				      */
    template <typename number2>
    void Tvmult (Vector<number2>&, const Vector<number2>&) const;

				     /**
				      * Execute backward block SOR
				      * preconditioning.
				      *
				      * Warning: this function
				      * performs normal @p vmult
				      * without adding. The reason for
				      * its existence is that
				      * BlockMatrixArray
				      * requires the adding version by
				      * default. On the other hand,
				      * adding requires an additional
				      * auxiliary vector, which is not
				      * desirable.
				      *
				      * @see vmult
				      */
    template <typename number2>
    void Tvmult_add (Vector<number2>&, const Vector<number2>&) const;

  private:
				     /**
				      * Actual implementation of the
				      * preconditioning algorithm.
				      *
				      * The parameter @p adding does
				      * not have any function, yet.
				      */
    template <typename number2>
    void do_vmult (Vector<number2>&,
		   const Vector<number2>&,
		   const bool adding) const;

				     /**
				      * Actual implementation of the
				      * preconditioning algorithm.
				      *
				      * The parameter @p adding does
				      * not have any function, yet.
				      */
    template <typename number2>
    void do_Tvmult (Vector<number2>&,
		    const Vector<number2>&,
		    const bool adding) const;

};


/**
 * Block SSOR preconditioning.
 *
 * The functions @p vmult and @p Tvmult execute a block-SSOR step,
 * based on the implementation in PreconditionBlockSOR.  This
 * class requires storage of the diagonal blocks and their inverses.
 *
 * See PreconditionBlock for requirements on the matrix.
 * @author Ralf Hartmann, Guido Kanschat, 1999, 2000
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class PreconditionBlockSSOR : public virtual Subscriptor,
			      private PreconditionBlockSOR<MATRIX, inverse_type>
{
  private:
				     /**
				      * Define number type of matrix.
				      */
    typedef typename MATRIX::value_type number;
    
  public:
				     /**
				      * Constructor.
				      */
    PreconditionBlockSSOR ();
				     /**
				      * Make initialization function
				      * publicly available.
				      */
    PreconditionBlockSOR<MATRIX,inverse_type>::initialize;
    
				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlockSOR<MATRIX,inverse_type>::clear;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlockSOR<MATRIX, inverse_type>::empty;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlockSOR<MATRIX, inverse_type>::el;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlockSOR<MATRIX,inverse_type>::set_same_diagonal;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlockSOR<MATRIX,inverse_type>::invert_diagblocks;

				     /**
				      * Execute block SSOR
				      * preconditioning.
				      *
				      * This function will
				      * automatically use the inverse
				      * matrices if they exist, if not
				      * then BlockSOR will waste much
				      * time inverting the diagonal
				      * block matrices in each
				      * preconditioning step.
				      */
    template <typename number2>
    void vmult (Vector<number2>&, const Vector<number2>&) const;

				     /**
				      * Same as vmult()
				      */
    template <typename number2>
    void Tvmult (Vector<number2>&, const Vector<number2>&) const;
};

/*@}*/
//----------------------------------------------------------------------//

template<class MATRIX, typename inverse_type>
inline
PreconditionBlock<MATRIX, inverse_type>::AdditionalData::
AdditionalData (const unsigned int block_size,
		const double relaxation,
		const bool invert_diagonal,
		const bool same_diagonal)
		:
		relaxation (relaxation),
		block_size(block_size),
		invert_diagonal(invert_diagonal),
		same_diagonal(same_diagonal)
{}


template<class MATRIX, typename inverse_type>
inline bool
PreconditionBlock<MATRIX, inverse_type>::empty () const
{
  if (A == 0)
    return true;
  return A->empty();
}


template<class MATRIX, typename inverse_type>
inline unsigned int
PreconditionBlock<MATRIX, inverse_type>::n_blocks () const
{
  return nblocks;
}


template<class MATRIX, typename inverse_type>
inline inverse_type
PreconditionBlock<MATRIX, inverse_type>::el (
  unsigned int i,
  unsigned int j) const
{
  const unsigned int bs = blocksize;
  const unsigned int nb = i/bs;
  
  const FullMatrix<inverse_type>& B = inverse(nb);

  const unsigned int ib = i % bs;
  const unsigned int jb = j % bs;

  if (jb + nb*bs != j)
    {
      return 0.;
    }
  
  return B(ib, jb);
}

//----------------------------------------------------------------------//

template<class MATRIX, typename inverse_type>
inline
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor::
Accessor (const PreconditionBlockJacobi<MATRIX, inverse_type>* matrix,
          const unsigned int row)
		:
		matrix(matrix),
		b_iterator(&matrix->inverse(0), 0, 0),
		b_end(&matrix->inverse(0), 0, 0)
{
  bs = matrix->block_size();
  a_block = row / bs;

				   // This is the end accessor, which
				   // does not hava a valid block.
  if (a_block == matrix->n_blocks())
    return;

  const unsigned int r = row % bs;

  b_iterator = matrix->inverse(a_block).begin(r);
  b_end = matrix->inverse(a_block).end();

  Assert (a_block < matrix->n_blocks(),
	  ExcIndexRange(a_block, 0, matrix->n_blocks()));
}


template<class MATRIX, typename inverse_type>
inline
unsigned int
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor::row() const
{
  Assert (a_block < matrix->n_blocks(),
	  ExcIteratorPastEnd());
  
  return bs * a_block + b_iterator->row();
}


template<class MATRIX, typename inverse_type>
inline
unsigned int
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor::column() const
{
  Assert (a_block < matrix->n_blocks(),
	  ExcIteratorPastEnd());
  
  return bs * a_block + b_iterator->column();
}


template<class MATRIX, typename inverse_type>
inline
inverse_type
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor::value() const
{
  Assert (a_block < matrix->n_blocks(),
	  ExcIteratorPastEnd());
  
  return b_iterator->value();
}


template<class MATRIX, typename inverse_type>
inline
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::
const_iterator(const PreconditionBlockJacobi<MATRIX, inverse_type> *matrix,
               const unsigned int row)
		:
		accessor(matrix, row)
{}


template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator &
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::operator++ ()
{
  Assert (*this != accessor.matrix->end(), ExcIteratorPastEnd());
  
  ++accessor.b_iterator;
  if (accessor.b_iterator == accessor.b_end)
    {
      ++accessor.a_block;
      
      if (accessor.a_block < accessor.matrix->n_blocks())
	{
	  accessor.b_iterator = accessor.matrix->inverse(accessor.a_block).begin();
	  accessor.b_end = accessor.matrix->inverse(accessor.a_block).end();
	}
    }
  return *this;
}


template<class MATRIX, typename inverse_type>
inline
const typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor &
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::operator* () const
{
  return accessor;
}


template<class MATRIX, typename inverse_type>
inline
const typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor *
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::operator-> () const
{
  return &accessor;
}


template<class MATRIX, typename inverse_type>
inline
bool
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::
operator == (const const_iterator& other) const
{
  if (accessor.a_block == accessor.matrix->n_blocks() &&
      accessor.a_block == other.accessor.a_block)
    return true;

  if (accessor.a_block != other.accessor.a_block)
    return false;
  
  return (accessor.row() == other.accessor.row() &&
          accessor.column() == other.accessor.column());
}


template<class MATRIX, typename inverse_type>
inline
bool
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::
operator != (const const_iterator& other) const
{
  return ! (*this == other);
}


template<class MATRIX, typename inverse_type>
inline
bool
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::
operator < (const const_iterator& other) const
{
  return (accessor.row() < other.accessor.row() ||
	  (accessor.row() == other.accessor.row() &&
           accessor.column() < other.accessor.column()));
}


template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator
PreconditionBlockJacobi<MATRIX, inverse_type>::begin () const
{
  return const_iterator(this, 0);
}


template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator
PreconditionBlockJacobi<MATRIX, inverse_type>::end () const
{
  return const_iterator(this, this->n_blocks() * this->block_size());
}


template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator
PreconditionBlockJacobi<MATRIX, inverse_type>::begin (
  const unsigned int r) const
{
  Assert (r < this->A->m(), ExcIndexRange(r, 0, this->A->m()));
  return const_iterator(this, r);
}



template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator
PreconditionBlockJacobi<MATRIX, inverse_type>::end (
  const unsigned int r) const
{
  Assert (r < this->A->m(), ExcIndexRange(r, 0, this->A->m()));
  return const_iterator(this, r+1);
}

#endif
