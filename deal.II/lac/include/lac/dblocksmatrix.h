/*----------------------------   dblocksmatrix.h     ---------------------------*/
/*      $Id$                 */
#ifndef __dblocksmatrix_H
#define __dblocksmatrix_H
/*----------------------------   dblocksmatrix.h     ---------------------------*/


#include <lac/dsmatrix.h>
#include <lac/dfmatrix.h>
#include <vector.h>

/**
 * Double precision block sparse matrix.
 * The block matrix assumes the matrix consisting of blocks on
 * the diagonal. These diagonal blocks and the elements below the
 * diagonal blocks are used in the #precondition_BlockSOR#.
 *
 * This block matrix structure is given e.g. for the DG method
 * for the transport equation and a downstream numbering.
 * If (as for this DG method) the matrix is empty above the
 * diagonal blocks BlockSOR is a direct solver.
 *
 * This first implementation of the BlockMatrix assumes the
 * matrix having blocks each of the same block size. Varying
 * block sizes within the matrix must still be implemented if needed.
 * @author Ralf Hartmann, 1999
 */
class dBlockSMatrix: public dSMatrix 
{
  public:
				     /**
				      * Constructor
				      */
    dBlockSMatrix();

				     /**
				      * Destructor
				      */
    virtual ~dBlockSMatrix();

				     /**
				      * Call #dSMatrix::reinit()# and
				      * delete the inverse matrices if existent.
				      */

    virtual void reinit();

				     /**
				      * Call #dSMatrix::reinit
				      * (const dSMatrixStruct &sparsity)# and
				      * delete the inverse matrices if existent.
				      */
    virtual void reinit (const dSMatrixStruct &sparsity);

				     /**
				      * Call #dSMatrix::clear# and 
				      * delete the inverse matrices if existent.
				      */
    virtual void clear ();

    				     /**
				      * Stores the inverse matrices of
				      * the diagonal blocks matrices
				      * in #inverse#. This costs some 
				      * additional memory (for DG
				      * methods about 1/3 of that used for
				      * the matrix) but it
				      * makes the preconditioning much faster.
				      */
    void invert_diagblocks();

				     /**
				      * Block SOR. Make sure that the right block size
				      * of the matrix is set by #set_block_size#
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
    void precondition_BlockSOR (dVector &dst, const dVector &src) const;
    
				     /**
				      * Set the right block size before calling
				      * #precondition_BlockSOR#.
				      * If block_size==1 BlockSOR is the same as SOR.
				      */
    void set_block_size (const unsigned int bsize);

				     /**
				      * Gives back the size of the blocks.
				      */
    unsigned int block_size() const;

				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongBlockSize,
		    int, int,
		    << "The blocksize " << arg1
		    << " and the size of the matrix " << arg2
		    << " do not match.");

    DeclException2 (ExcWrongInverses,
		    int, int,
		    << "There are " << arg1
		    << " inverse matrices but " << arg2
		    << " cells.");

				     /**
				      * Exception
				      */
    DeclException0 (ExcInverseMatricesDoNotExist);

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
    DeclException0 (ExcInternalError);
    
  private:
				     /**
				      * size of the blocks.
				      */
    unsigned int blocksize;

				     /**
				      * stores the inverse matrices of
				      * the diagonal blocks matrices
				      */
    vector<dFMatrix> inverse;
};



/*----------------------------   dblocksmatrix.h     ---------------------------*/
/* end of #ifndef __dblocksmatrix_H */
#endif
/*----------------------------   dblocksmatrix.h     ---------------------------*/
