//----------------------------  precondition_block.h  ---------------------------
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
//----------------------------  precondition_block.h  ---------------------------
#ifndef __deal2__precondition_block_h
#define __deal2__precondition_block_h


#include <lac/forward_declarations.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>

#include <vector>


/**
 * Base class for #PreconditionBlockJacobi#, #PreconditionBlockSOR#, ...
 * This class assumes the #SparseMatrix<number># consisting of invertible blocks 
 * of #blocksize# on the diagonal and provides the inversion of the diagonal blocks
 * of the matrix. NOT only block diagonal matrices are allowed but all
 * matrices of arbitrary structure with the minimal property of having
 * invertible blocks on the diagonal!
 *
 * This block matrix structure is given e.g. for the DG method
 * for the transport equation. For a downstream numbering the matrices
 * even have got a block lower left matrix structure, i.e. the matrices
 * are empty above the diagonal blocks.
 *
 * For all matrices that are empty above and below the diagonal
 * blocks (i.e. for all block diagonal matrices) the #BlockJacobi# preconditioner
 * is a direct solver. For all matrices that are empty only above the diagonal blocks
 * (e.g. the matrices one gets by the DG method with downstream numbering) the
 * #BlockSOR# is a direct solver.
 * 
 * This first implementation of the #PreconditionBlock# assumes the
 * matrix has blocks each of the same block size. Varying
 * block sizes within the matrix must still be implemented if needed.
 *
 * The first template parameter denotes the type of number representation in
 * the sparse matrix, the second denotes the type of number representation in
 * which the inverted diagonal block matrices are stored within this class
 * by #invert_diagblocks()#. If you don't want to use the block inversion as
 * an exact solver, but rather as a preconditioner, you may probably want to
 * store the inverted blocks with less accuracy than the original matrix;
 * for example, #number==double, inverse_type=float# might be a viable choice.
 *
 *
 * \section{On template instantiations}
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
 * @author Ralf Hartmann, 1999
 */
template<typename number, typename inverse_type = number>
class PreconditionBlock: public Subscriptor
{
  public:
				     /**
				      * Constructor.
				      */
    PreconditionBlock();
    
				     /**
				      * Destructor.
				      */
    virtual ~PreconditionBlock();

				     /**
				      * Deletes the inverse diagonal block
				      * matrices if existent, sets the
				      * blocksize to 0, hence leaves the
				      * class in the state that it had 
				      * directly after
				      * calling the constructor.
				      */
    virtual void clear();

				     /**
				      * Takes the matrix that should be used
				      * for the preconditioning. A reference
				      * to it is stored within this class,
				      * but ownership of the matrix remains
				      * with the caller of this function.
				      */
    void use_matrix(const SparseMatrix<number> &M);
    
				     /**
				      * Set the right block size before calling
				      * #invert_diagblocks# or calling the
				      * #operator ()# function
				      * of #PreconditionBlockSOR#
				      * or #PreconditionBlockJacobi#.
				      * If #block_size==1# BlockSOR or 
				      * BlockJacobi are equal to the
				      * standard SOR or Jacobi 
				      * preconditioner.
				      */
    void set_block_size (const unsigned int bsize);

    				     /**
				      * Stores the inverse of
				      * the diagonal blocks
				      * in #inverse#. This costs some 
				      * additional memory - for DG
				      * methods about 1/3 (for double inverses) 
				      * or 1/6 (for float inverses) of that
				      * used for the matrix - but it
				      * makes the preconditioning much faster.
				      *
				      * It is not allowed to call this function
				      * twice (will produce an error) before
				      * a call of #clear(..)#
				      * because at the second time there already
				      * exist the inverse matrices.
				      *
				      * After this function is called, the
				      * lock on the matrix given through the
				      * #use_matrix# function is released,
				      * i.e. you may overwrite of delete it.
				      * You may want to do this in case
				      * you use this matrix to precondition
				      * another matrix.
				      */
    void invert_diagblocks();

				     /**
				      * Return the size of the blocks.
				      */
    unsigned int block_size () const;

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
    DeclException0 (ExcBlockSizeNotSet);

				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotSquare);

				     /**
				      * Exception
				      */
    DeclException0 (ExcNoMatrixGivenToUse);
   
  protected:
				     /**
				      * Size of the blocks. Each diagonal
				      * block is assumed to be of the
				      * same size.
				      */
    unsigned int blocksize;
    
				     /**
				      * Storage of the inverse matrices of
				      * the diagonal blocks matrices as
				      * #FullMatrix<inverse_type># matrices.
				      * For BlockSOR as preconditioning
				      * using #inverse_type=float# saves memory
				      * in comparison with #inverse_type=double#.
				      */
    vector<FullMatrix<inverse_type> > inverse;

				     /**
				      * Pointer to the matrix. Make sure that
				      * the matrix exists as long as this class
				      * needs it, i.e. until calling #invert_diagblocks#,
				      * or (if the inverse matrices should not be
				      * stored) until the last call of the 
				      * preconditoining #operator()# function of the
				      * derived classes.
				      */
    SmartPointer<const SparseMatrix<number> > A;
};


/**
 * Block SOR preconditioning.
 *
 * The diagonal blocks and the elements 
 * (of arbitray structure) below the diagonal blocks are used
 * in the #operator ()# function of this class.
 */
template<typename number,
         typename inverse_type = number>
class PreconditionBlockSOR : public PreconditionBlock<number,inverse_type>
{
  public:
				     /**
				      * Constructor. Takes the damping
				      * Parameter as 
				      */
    PreconditionBlockSOR (const number omega=1.);
    
				     /**
				      * Destructor.
				      */
    virtual ~PreconditionBlockSOR();

				     /**
				      * Set the relaxation parameter.
				      */
    void set_omega(number omega);
    
				     /**
				      * Execute block SOR preconditioning.
				      * Make sure that the right
				      * block size
				      * of the matrix is set by
				      * #set_block_size#
				      * before calling this function.
				      *
				      * This function will automatically use the
				      * inverse matrices if they exist, if not
				      * then BlockSOR will waste much time
				      * inverting the diagonal block
				      * matrices in each preconditioning step.
				      *
				      * For matrices which are
				      * empty above the diagonal blocks
				      * BlockSOR is a direct solver.
				      */
    template <typename number2>
    void operator() (Vector<number2>&, const Vector<number2>&) const;

  private:
				     /**
				      * Damping parameter.
				      */
    number omega;
};


#endif

