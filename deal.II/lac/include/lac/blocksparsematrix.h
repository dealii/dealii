/*----------------------------   blocksparsematrix.h     ---------------------*/
/*      $Id$                 */
/*           Ralf Hartmann, University of Heidelberg                          */
#ifndef __blocksparsematrix_H
#define __blocksparsematrix_H
/*----------------------------   blocksparsematrix.h     ---------------------*/


#include <lac/sparsematrix.h>
#include <lac/forward-declarations.h>

#include <vector>




/**
 * Block sparse matrix.
 * The block matrix assumes the matrix consisting of blocks on
 * the diagonal. These diagonal blocks and the elements 
 * (of arbitray structure) below the
 * diagonal blocks are used in the #precondition_BlockSOR# function.
 *
 * This block matrix structure is given e.g. for the DG method
 * for the transport equation and a downstream numbering.
 * If (as for this DG method) the matrix is empty above the
 * diagonal blocks then BlockSOR is a direct solver.
 *
 * This first implementation of the BlockMatrix assumes the
 * matrix having blocks each of the same block size. Varying
 * block sizes within the matrix must still be implemented if needed.
 *
 * The first template parameter denotes the type of number representation in
 * the sparse matrix, the second denotes the type of number representation in
 * which the inverse diagonal block
 * matrices are stored by #invert_diagblocks()#.
 * 
 * @author Ralf Hartmann, 1999
 */
template <typename number, typename inverse_type>
class BlockSparseMatrix: public SparseMatrix<number> 
{
  public:
				     /**
				      * Constructor
				      */
    BlockSparseMatrix();

				     /**
				      * Destructor
				      */
    virtual ~BlockSparseMatrix();

				     /**
				      * Call #dSMatrix::reinit()# and
				      * delete the inverse diagonal block matrices if existent.
				      */
    virtual void reinit();

				     /**
				      * Call #dSMatrix::reinit
				      * (const dSMatrixStruct &sparsity)# and
				      * delete the inverse diagonal block matrices if existent.
				      */
    virtual void reinit (const SparseMatrixStruct &sparsity);

				     /**
				      * Call #dSMatrix::clear# and 
				      * delete the inverse diagonal block matrices if existent.
				      */
    virtual void clear ();

    				     /**
				      * Stores the inverse of
				      * the diagonal blocks matrices
				      * in #inverse#. This costs some 
				      * additional memory - for DG
				      * methods about 1/3 (for double inverses) 
				      * or 1/6 (for float inverses) of that
				      * used for the matrix - but it
				      * makes the preconditioning much faster.
				      *
				      * It is not allowed to call this function
				      * twice (will produce an error) before
				      * a call of #reinit(..)#
				      * because at the second time there already
				      * exist the inverse matrices.
				      */
    void invert_diagblocks();

				     /**
				      * Block SOR. Make sure that the right
				      * block size
				      * of the matrix is set by
				      * #set_block_size#
				      * before calling this function.
				      *
				      * BlockSOR will automatically use the
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
    void precondition_BlockSOR (Vector<number2> &dst, const Vector<number2> &src,
				const number om = 1.) const;
    
				     /**
				      * Set the right block size before calling
				      * #precondition_BlockSOR#.
				      * If block_size==1 BlockSOR is the
				      * same as SOR.
				      */
    void set_block_size (const unsigned int bsize);

				     /**
				      * Gives back the size of the blocks.
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
   
  private:
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
};


/*----------------------------   blocksparsematrix.h     ---------------------------*/
/* end of #ifndef __blocksparsematrix_H */
#endif
/*----------------------------   blocksparsematrix.h     ---------------------------*/
