//----------------------------  precondition_block.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  precondition_block.h  ---------------------------
#ifndef __deal2__precondition_block_h
#define __deal2__precondition_block_h


#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>

#include <vector>

template <typename number> class FullMatrix;
template <typename number> class SparseMatrix;
template <typename number> class Vector;


/**
 * Base class for @p{PreconditionBlockJacobi},
 * @p{PreconditionBlockSOR}, ...  This class assumes the
 * @p{SparseMatrix<number>} consisting of invertible blocks of
 * @p{blocksize} on the diagonal and provides the inversion of the
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
 * (i.e. for all block diagonal matrices) the @p{BlockJacobi}
 * preconditioner is a direct solver. For all matrices that are empty
 * only above the diagonal blocks (e.g. the matrices one gets by the
 * DG method with downstream numbering) @p{BlockSOR} is a direct
 * solver.
 * 
 * This first implementation of the @p{PreconditionBlock} assumes the
 * matrix has blocks each of the same block size. Varying block sizes
 * within the matrix must still be implemented if needed.
 *
 * The first template parameter denotes the type of number
 * representation in the sparse matrix, the second denotes the type of
 * number representation in which the inverted diagonal block matrices
 * are stored within this class by @p{invert_diagblocks()}. If you
 * don't want to use the block inversion as an exact solver, but
 * rather as a preconditioner, you may probably want to store the
 * inverted blocks with less accuracy than the original matrix; for
 * example, @p{number==double, inverse_type=float} might be a viable
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
template<typename number, typename inverse_type = number>
class PreconditionBlock
{
  public:
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
    void initialize (const SparseMatrix<number>& A,
		     const unsigned int block_size,
		     const double relaxation = 1.);

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
				      * Stores the inverse of the
				      * diagonal blocks in
				      * @p{inverse}. This costs some
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
				      * @p{clear(...)}  because at the
				      * second time there already
				      * exist the inverse matrices.
				      *
				      * After this function is called,
				      * the lock on the matrix given
				      * through the @p{use_matrix}
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
   
  protected:
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
				      * @p{invert_diagblocks}, or (if
				      * the inverse matrices should
				      * not be stored) until the last
				      * call of the preconditoining
				      * @p{vmult} function of the
				      * derived classes.
				      */
    SmartPointer<const SparseMatrix<number> > A;
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

  private:
				     /**
				      * Storage of the inverse
				      * matrices of the diagonal
				      * blocks matrices as
				      * @p{FullMatrix<inverse_type>}
				      * matrices. Using
				      * @p{inverse_type=float} saves
				      * memory in comparison with
				      * @p{inverse_type=double}.
				      */
    typename std::vector<FullMatrix<inverse_type> > var_inverse;

				     /**
				      * Storage of the original diagonal blocks.
				      * These are only filled if @p{store_diagonals}
				      * is @p{true}.
				      *
				      * Used by the blocked SSOR method.
				      */
    typename std::vector<FullMatrix<inverse_type> > var_diagonal;
				      
				     /**
				      * Flag for diagonal compression.
				      * @see set_same_diagonal()
				      */
    bool same_diagonal;
};



/**
 * Block Jacobi preconditioning.
 * See @ref{PreconditionBlock} for requirements on the matrix.
 * @author Ralf Hartmann, Guido Kanschat, 1999, 2000
 */
template<typename number, typename inverse_type = number>
class PreconditionBlockJacobi : public Subscriptor,
				private PreconditionBlock<number,inverse_type>
{
  public:
				     /**
				      * Make initialization function
				      * publicly available.
				      */
    PreconditionBlock<number,inverse_type>::initialize;
    
				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<number,inverse_type>::clear;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<number,inverse_type>::set_same_diagonal;

				     /**
				      * Make function public.
				      */
    PreconditionBlock<number,inverse_type>::invert_diagblocks;
    
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
				      * Same as @p{vmult}, since Jacobi is symmetric.
				      */
    template <typename number2>
    void Tvmult (Vector<number2>&, const Vector<number2>&) const;
};



/**
 * Block SOR preconditioning.
 *
 * The functions @p{vmult} and @p{Tvmult} execute a (transposed)
 * block-SOR step, based on the blocks in @ref{PreconditionBlock}. The
 * elements outside the diagonal blocks may be distributed
 * arbitrarily.
 *
 * See @ref{PreconditionBlock} for requirements on the matrix.
 * @author Ralf Hartmann, Guido Kanschat, 1999, 2000
 */
template<typename number, typename inverse_type = number>
class PreconditionBlockSOR : public  Subscriptor,
			     protected PreconditionBlock<number,inverse_type>
{
  public:
				     /**
				      * Make initialization function
				      * publicly available.
				      */
    PreconditionBlock<number,inverse_type>::initialize;
    
				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<number,inverse_type>::clear;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<number,inverse_type>::set_same_diagonal;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<number,inverse_type>::invert_diagblocks;

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
				      * Transpose of @ref{vmult}
				      */
    template <typename number2>
    void Tvmult (Vector<number2>&, const Vector<number2>&) const;
};


/**
 * Block SSOR preconditioning.
 *
 * The functions @p{vmult} and @p{Tvmult} execute a block-SSOR step,
 * based on the implementation in @ref{PreconditionBlockSOR}.  This
 * class requires storage of the diagonal blocks and their inverses.
 *
 * See @ref{PreconditionBlock} for requirements on the matrix.
 * @author Ralf Hartmann, Guido Kanschat, 1999, 2000
 */
template<typename number, typename inverse_type = number>
class PreconditionBlockSSOR : private PreconditionBlockSOR<number,inverse_type>
{
  public:
				     /**
				      * Constructor.
				      */
    PreconditionBlockSSOR ();
				     /**
				      * Make initialization function
				      * publicly available.
				      */
    PreconditionBlock<number,inverse_type>::initialize;
    
				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<number,inverse_type>::clear;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<number,inverse_type>::set_same_diagonal;

				     /**
				      * Make function of base class public again.
				      */
    PreconditionBlock<number,inverse_type>::invert_diagblocks;

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
				      * Same as @ref{vmult}
				      */
    template <typename number2>
    void Tvmult (Vector<number2>&, const Vector<number2>&) const;
};


#endif

