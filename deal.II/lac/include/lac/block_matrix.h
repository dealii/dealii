//-------------------------------------------------------------------
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
//-------------------------------------------------------------------
#ifndef __deal2__block_matrix_h
#define __deal2__block_matrix_h


#include <base/exceptions.h>
#include <base/smartpointer.h>
#include <lac/block_vector.h>

/**
 * A matrix with several copies of the same block on the diagonal.
 *
 * This matrix implements an @p{m} by @p{m} block matrix. Each
 * diagonal block consists of the same (non-block) matrix, while
 * off-diagonal blocks are void.
 *
 * One special application is a one by one block matrix, allowing to
 * apply the @p{vmult} of the original matrix (or preconditioner) to a
 * block vector.
 *
 * @author Guido Kanschat, 2000
 */
template <class MATRIX>
class BlockDiagonalMatrix : public Subscriptor
{
public:
				   /**
				    * Constructor for an @p{n_blocks}
				    * by @p{n_blocks} matrix with
				    * diagonal blocks @p{M}.
				    */
  BlockDiagonalMatrix (const MATRIX& M,
		       unsigned int n_blocks);

				   /**
				    * Matrix-vector-multiplication.
				    */
  template <typename number1, typename number2>
  void vmult (BlockVector<number1>& dst,
	      const BlockVector<number2>& src) const;
  
				   /**
				    * Transposed matrix-vector-multiplication.
				    */
  template <typename number1, typename number2>
  void Tvmult (BlockVector<number1>& dst,
	       const BlockVector<number2>& src) const;
private:
				   /**
				    * Number of blocks.
				    */
  unsigned int num_blocks;
				   /**
				    * Diagonal entry.
				    */
  SmartPointer<const MATRIX> matrix;
};


//----------------------------------------------------------------------//

template <class MATRIX>
BlockDiagonalMatrix<MATRIX>::BlockDiagonalMatrix (const MATRIX& M,
						  unsigned int num_blocks)
  : num_blocks (num_blocks),
    matrix(&M)
{}


template <class MATRIX>
template <typename number1, typename number2>
void
BlockDiagonalMatrix<MATRIX>::vmult (BlockVector<number1>& dst,
				     const BlockVector<number2>& src) const
{
  Assert (dst.n_blocks()==num_blocks,
	  ExcDimensionMismatch(dst.n_blocks(),num_blocks));
  Assert (src.n_blocks()==num_blocks,
	  ExcDimensionMismatch(src.n_blocks(),num_blocks));

  for (unsigned int i=0;i<num_blocks;++i)
    matrix->vmult (dst.block(i), src.block(i));
}


template <class MATRIX>
template <typename number1, typename number2>
void
BlockDiagonalMatrix<MATRIX>::Tvmult (BlockVector<number1>& dst,
				     const BlockVector<number2>& src) const
{
  Assert (dst.n_blocks()==num_blocks,
	  ExcDimensionMismatch(dst.n_blocks(),num_blocks));
  Assert (src.n_blocks()==num_blocks,
	  ExcDimensionMismatch(src.n_blocks(),num_blocks));

  for (unsigned int i=0;i<num_blocks;++i)
    matrix->Tvmult (dst.block(i), src.block(i));
}


#endif

