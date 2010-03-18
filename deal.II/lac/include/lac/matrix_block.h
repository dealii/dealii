//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2009, 2010 by the deal.II authors
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
#include <base/named_data.h>

#include <boost/shared_ptr.hpp>

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
 * @code
 * std::vector<MatrixBlock<SparseMatrix<double> > > blocks;
 *
 * ...
 *
 * BlockMatrixArray matrix (n_blocks, n_blocks);
 *
 * for (unsigned int i=0;i<blocks.size;++i)
 *   matrix.enter(blocks[i].row, blocks[i].column, blocks[i].matrix);
 * @endcode
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


/**
 * A vector of MatrixBlock, which is implemented using shared
 * pointers, in order to allow for copying and rearranging. Each
 * matrix block can be identified by name.
 *
 * @author Baerbel Janssen, Guido Kanschat, 2010
 */
template <class MATRIX>
class MatrixBlockVector : public NamedData<boost::shared_ptr<MatrixBlock<MATRIX> > >
{
  public:
				     /**
				      * Add a new matrix block at the
				      * position <tt(row,column)</tt>
				      * in the block system.
				      */
    void add(unsigned int row, unsigned int column, const std::string& name);
				     /**
				      * Access a constant reference to
				      * the block at position <i>i</i>.
				      */
    const MatrixBlock<MATRIX>& block(unsigned int i) const;
				     /**
				      * Access the matrix at position
				      * <i>i</i> for read and write
				      * access.
				      */
    MATRIX& matrix(unsigned int i);
    
    
};


//----------------------------------------------------------------------//

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

//----------------------------------------------------------------------//

template <class MATRIX>
inline void
MatrixBlockVector<MATRIX>::add(unsigned int row, unsigned int column, const std::string& name)
{
  boost::shared_ptr<MatrixBlock<MATRIX> > p(new MatrixBlock<MATRIX>(row, column));
  NamedData<boost::shared_ptr<MatrixBlock<MATRIX> > >::add(p, name);
}


template <class MATRIX>
inline const MatrixBlock<MATRIX>&
MatrixBlockVector<MATRIX>::block(unsigned int i) const
{
  return *this->read(i);
}


template <class MATRIX>
inline MATRIX&
MatrixBlockVector<MATRIX>::matrix(unsigned int i)
{
  return (*this)(i)->matrix;
}



DEAL_II_NAMESPACE_CLOSE

#endif
