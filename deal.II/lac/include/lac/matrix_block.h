//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__matrix_block_h
#define __deal2__matrix_block_h

#include <base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A wrapper around a matrix object, storing the coordinates in a
 * block matrix as well.
 *
 * @note They are the position in the global matrix, where we don't use
 * BlockSparseMatrix. Matrices may have the same coordinates for the
 * following reason: A preconditioner for the Oseen system can be
 * built as a block system, where the pressure block is of the form
 * M^-1 F A^-1 with m the mass matrix (u,v), A the Laplacian and F the
 * advection diffusion operator. So, I need to build three matrices
 * for a single block in my system + in some cases a pressure
 * stabilization.
 *
 * MatrixBlock comes handy when using BlockMatrixArray. Once the
 * MatrixBlock has been properly initalized and filled, it can be used
 * as:
 * @begin{code}
 * std::vector<MatrixBlock<SparseMatrix<double> > > blocks;
 *
 * ...
 *
 * BlockMatrixArray matrix (n_blocks, n_blocks);
 *
 * for (unsigned int i=0;i<blocks.size;++i)
 *   matrix.enter(blocks[i].row, blocks[i].column, blocks[i].matrix);
 * @end{code}
 *
 *
 * @author Guido Kanschat, 2006
 */
template <class MATRIX>
struct MatrixBlock
{
				     /**
				      * Constructor rendering an
				      * uninitialized object.
				      */
    MatrixBlock();
    
				     /**
				      * Copy constructor.
				      */
    MatrixBlock(const MatrixBlock<MATRIX>& M);
    
				     /**
				      * Constructor setting block
				      * coordinates, but not
				      * initializing the matrix.
				      */
    
    MatrixBlock(unsigned int i, unsigned int j);
    
				     /**
				      * Row coordinate.  This is the
				      * position of the data member
				      * matrix on the global matrix.
				      */
    unsigned int row;
				     /**
				      * Column coordinate.  This is
				      * the position of the data
				      * member matrix on the global
				      * matrix.
				      */
    unsigned int column;
    
				     /**
				      * The matrix itself
				      */
    MATRIX matrix;
};


template <class MATRIX>
MatrixBlock<MATRIX>::MatrixBlock()
		:
		row(deal_II_numbers::invalid_unsigned_int),
		column(deal_II_numbers::invalid_unsigned_int)
{}


template <class MATRIX>
MatrixBlock<MATRIX>::MatrixBlock(const MatrixBlock<MATRIX>& M)
		:
		row(M.row),
		column(M.column),
		matrix(M.matrix)
{}


template <class MATRIX>
MatrixBlock<MATRIX>::MatrixBlock(unsigned int i, unsigned int j)
		:
		row(i), column(j)
{}


DEAL_II_NAMESPACE_CLOSE

#endif
