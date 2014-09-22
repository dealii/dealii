// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__block_matrix_h
#define __deal2__block_matrix_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/block_vector.h>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Matrix2
 *@{
 */

/**
 * A matrix with several copies of the same block on the diagonal.
 *
 * This matrix implements an @p m by @p m block matrix. Each
 * diagonal block consists of the same (non-block) matrix, while
 * off-diagonal blocks are void.
 *
 * One special application is a one by one block matrix, allowing to
 * apply the @p vmult of the original matrix (or preconditioner) to a
 * block vector.
 *
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Guido Kanschat, 2000
 */
template <class MATRIX>
class BlockDiagonalMatrix : public Subscriptor
{
public:
  /**
   * Constructor for an @p n_blocks
   * by @p n_blocks matrix with
   * diagonal blocks @p M.
   */
  BlockDiagonalMatrix (const MATRIX       &M,
                       const unsigned int  n_blocks);

  /**
   * Matrix-vector-multiplication.
   */
  template <typename number1, typename number2>
  void vmult (BlockVector<number1> &dst,
              const BlockVector<number2> &src) const;

  /**
   * Transposed matrix-vector-multiplication.
   */
  template <typename number1, typename number2>
  void Tvmult (BlockVector<number1> &dst,
               const BlockVector<number2> &src) const;
private:
  /**
   * Number of blocks.
   */
  unsigned int num_blocks;

  /**
   * Diagonal entry.
   */
  SmartPointer<const MATRIX,BlockDiagonalMatrix<MATRIX> > matrix;
};

/*@}*/
//---------------------------------------------------------------------------

template <class MATRIX>
BlockDiagonalMatrix<MATRIX>::BlockDiagonalMatrix (const MATRIX &M,
                                                  const unsigned int num_blocks)
  :
  num_blocks (num_blocks),
  matrix(&M)
{}


template <class MATRIX>
template <typename number1, typename number2>
void
BlockDiagonalMatrix<MATRIX>::vmult (BlockVector<number1> &dst,
                                    const BlockVector<number2> &src) const
{
  Assert (dst.n_blocks()==num_blocks,
          ExcDimensionMismatch(dst.n_blocks(),num_blocks));
  Assert (src.n_blocks()==num_blocks,
          ExcDimensionMismatch(src.n_blocks(),num_blocks));

  for (unsigned int i=0; i<num_blocks; ++i)
    matrix->vmult (dst.block(i), src.block(i));
}


template <class MATRIX>
template <typename number1, typename number2>
void
BlockDiagonalMatrix<MATRIX>::Tvmult (BlockVector<number1> &dst,
                                     const BlockVector<number2> &src) const
{
  Assert (dst.n_blocks()==num_blocks,
          ExcDimensionMismatch(dst.n_blocks(),num_blocks));
  Assert (src.n_blocks()==num_blocks,
          ExcDimensionMismatch(src.n_blocks(),num_blocks));

  for (unsigned int i=0; i<num_blocks; ++i)
    matrix->Tvmult (dst.block(i), src.block(i));
}


DEAL_II_NAMESPACE_CLOSE

#endif
