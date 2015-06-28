// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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


#include <deal.II/lac/block_matrix_array.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>

DEAL_II_NAMESPACE_OPEN


template <typename number, typename BLOCK_VECTOR>
BlockMatrixArray<number,BLOCK_VECTOR>::Entry::Entry (const Entry &e)
  :
  row(e.row),
  col(e.col),
  prefix(e.prefix),
  transpose(e.transpose),
  matrix(e.matrix)
{
  Entry &e2 = const_cast<Entry &>(e);
  e2.matrix = 0;
}



template <typename number, typename BLOCK_VECTOR>
BlockMatrixArray<number,BLOCK_VECTOR>::Entry::~Entry ()
{
  if (matrix)
    delete matrix;
}


template <typename number, typename BLOCK_VECTOR>
BlockMatrixArray<number,BLOCK_VECTOR>::BlockMatrixArray ()
  : block_rows (0),
    block_cols (0)
{}



template <typename number, typename BLOCK_VECTOR>
BlockMatrixArray<number,BLOCK_VECTOR>::BlockMatrixArray (
  const unsigned int n_block_rows,
  const unsigned int n_block_cols)
  : block_rows (n_block_rows),
    block_cols (n_block_cols)
{}


template <typename number, typename BLOCK_VECTOR>
void
BlockMatrixArray<number,BLOCK_VECTOR>::initialize (
  const unsigned int n_block_rows,
  const unsigned int n_block_cols)
{
  block_rows = n_block_rows;
  block_cols = n_block_cols;
}



template <typename number, typename BLOCK_VECTOR>
void
BlockMatrixArray<number,BLOCK_VECTOR>::reinit (
  const unsigned int n_block_rows,
  const unsigned int n_block_cols)
{
  clear();
  block_rows = n_block_rows;
  block_cols = n_block_cols;
}



template <typename number, typename BLOCK_VECTOR>
void
BlockMatrixArray<number,BLOCK_VECTOR>::clear ()
{
  entries.clear();
}


template <typename number, typename BLOCK_VECTOR>
void
BlockMatrixArray<number,BLOCK_VECTOR>::vmult_add (BLOCK_VECTOR &dst,
                                                  const BLOCK_VECTOR &src) const
{
  GrowingVectorMemory<typename BLOCK_VECTOR::BlockType > mem;
  Assert (dst.n_blocks() == block_rows,
          ExcDimensionMismatch(dst.n_blocks(), block_rows));
  Assert (src.n_blocks() == block_cols,
          ExcDimensionMismatch(src.n_blocks(), block_cols));

  typename VectorMemory<typename BLOCK_VECTOR::BlockType >::Pointer p_aux(mem);
  typename BLOCK_VECTOR::BlockType &aux = *p_aux;

  typename std::vector<Entry>::const_iterator m = entries.begin();
  typename std::vector<Entry>::const_iterator end = entries.end();

  for (; m != end ; ++m)
    {
      aux.reinit(dst.block(m->row));
      if (m->transpose)
        m->matrix->Tvmult(aux, src.block(m->col));
      else
        m->matrix->vmult(aux, src.block(m->col));
      dst.block(m->row).add (m->prefix, aux);
    }
}




template <typename number, typename BLOCK_VECTOR>
void
BlockMatrixArray<number,BLOCK_VECTOR>::vmult (BLOCK_VECTOR &dst,
                                              const BLOCK_VECTOR &src) const
{
  dst = 0.;
  vmult_add (dst, src);
}




template <typename number, typename BLOCK_VECTOR>
void
BlockMatrixArray<number,BLOCK_VECTOR>::Tvmult_add (BLOCK_VECTOR &dst,
                                                   const BLOCK_VECTOR &src) const
{
  GrowingVectorMemory<typename BLOCK_VECTOR::BlockType > mem;
  Assert (dst.n_blocks() == block_cols,
          ExcDimensionMismatch(dst.n_blocks(), block_cols));
  Assert (src.n_blocks() == block_rows,
          ExcDimensionMismatch(src.n_blocks(), block_rows));

  typename std::vector<Entry>::const_iterator m = entries.begin();
  typename std::vector<Entry>::const_iterator end = entries.end();

  typename VectorMemory<typename BLOCK_VECTOR::BlockType >::Pointer p_aux(mem);
  typename BLOCK_VECTOR::BlockType &aux = *p_aux;

  for (; m != end ; ++m)
    {
      aux.reinit(dst.block(m->col));
      if (m->transpose)
        m->matrix->vmult(aux, src.block(m->row));
      else
        m->matrix->Tvmult(aux, src.block(m->row));
      dst.block(m->col).add (m->prefix, aux);
    }
}



template <typename number, typename BLOCK_VECTOR>
void
BlockMatrixArray<number,BLOCK_VECTOR>::Tvmult (BLOCK_VECTOR &dst,
                                               const BLOCK_VECTOR &src) const
{
  dst = 0.;
  Tvmult_add (dst, src);
}




template <typename number, typename BLOCK_VECTOR>
number
BlockMatrixArray<number,BLOCK_VECTOR>::matrix_scalar_product (
  const BLOCK_VECTOR &u,
  const BLOCK_VECTOR &v) const
{
  GrowingVectorMemory<typename BLOCK_VECTOR::BlockType > mem;
  Assert (u.n_blocks() == block_rows,
          ExcDimensionMismatch(u.n_blocks(), block_rows));
  Assert (v.n_blocks() == block_cols,
          ExcDimensionMismatch(v.n_blocks(), block_cols));

  typename VectorMemory<typename BLOCK_VECTOR::BlockType >::Pointer p_aux(mem);
  typename BLOCK_VECTOR::BlockType &aux = *p_aux;

  typename std::vector<Entry>::const_iterator m;
  typename std::vector<Entry>::const_iterator end = entries.end();

  number result = 0.;

  for (unsigned int i=0; i<block_rows; ++i)
    {
      aux.reinit(u.block(i));
      for (m = entries.begin(); m != end ; ++m)
        {
          if (m->row != i)
            continue;
          if (m->transpose)
            m->matrix->Tvmult_add(aux, v.block(m->col));
          else
            m->matrix->vmult(aux, v.block(m->col));
        }
      result += u.block(i)*aux;
    }

  return result;
}



template <typename number, typename BLOCK_VECTOR>
number
BlockMatrixArray<number,BLOCK_VECTOR>::matrix_norm_square (
  const BLOCK_VECTOR &u) const
{
  return matrix_scalar_product(u,u);
}



template <typename number, typename BLOCK_VECTOR>
unsigned int
BlockMatrixArray<number,BLOCK_VECTOR>::n_block_rows () const
{
  return block_rows;
}



template <typename number, typename BLOCK_VECTOR>
unsigned int
BlockMatrixArray<number,BLOCK_VECTOR>::n_block_cols () const
{
  return block_cols;
}



//---------------------------------------------------------------------------

template <typename number, typename BLOCK_VECTOR>
BlockTrianglePrecondition<number,BLOCK_VECTOR>::BlockTrianglePrecondition()
  : BlockMatrixArray<number,BLOCK_VECTOR> (),
    backward(false)
{}


template <typename number, typename BLOCK_VECTOR>
BlockTrianglePrecondition<number,BLOCK_VECTOR>::BlockTrianglePrecondition(
  const unsigned int block_rows)
  :
  BlockMatrixArray<number,BLOCK_VECTOR> (block_rows, block_rows),
  backward(false)
{}


template <typename number, typename BLOCK_VECTOR>
void
BlockTrianglePrecondition<number,BLOCK_VECTOR>::reinit (
  const unsigned int n)
{
  BlockMatrixArray<number,BLOCK_VECTOR>::reinit(n,n);
}


template <typename number, typename BLOCK_VECTOR>
void
BlockTrianglePrecondition<number,BLOCK_VECTOR>::do_row (
  BLOCK_VECTOR &dst,
  size_type row_num) const
{
  GrowingVectorMemory<typename BLOCK_VECTOR::BlockType > mem;
  typename std::vector<typename BlockMatrixArray<number,BLOCK_VECTOR>::Entry>::const_iterator
  m = this->entries.begin();
  typename std::vector<typename BlockMatrixArray<number,BLOCK_VECTOR>::Entry>::const_iterator
  end = this->entries.end();
  std::vector<typename std::vector<typename BlockMatrixArray<number,BLOCK_VECTOR>::Entry>::const_iterator>
  diagonals;

  typename VectorMemory<typename BLOCK_VECTOR::BlockType >::Pointer p_aux(mem);
  typename BLOCK_VECTOR::BlockType &aux = *p_aux;

  aux.reinit(dst.block(row_num), true);

  // Loop over all entries, since
  // they are not ordered by rows.
  for (; m != end ; ++m)
    {
      const size_type i=m->row;
      // Ignore everything not in
      // this row
      if (i != row_num)
        continue;
      const size_type j=m->col;
      // Only use the lower (upper)
      // triangle for forward
      // (backward) substitution
      if (((j > i) && !backward) || ((j < i) && backward))
        continue;
      if (j == i)
        {
          diagonals.push_back(m);
        }
      else
        {
          if (m->transpose)
            m->matrix->Tvmult(aux, dst.block(j));
          else
            m->matrix->vmult(aux, dst.block(j));
          dst.block(i).add (-1 * m->prefix, aux);
        }
    }
  Assert (diagonals.size() != 0, ExcNoDiagonal(row_num));

  // Inverting the diagonal block is
  // simple, if there is only one
  // matrix
  if (diagonals.size() == 1)
    {
      if (diagonals[0]->transpose)
        diagonals[0]->matrix->Tvmult(aux, dst.block(row_num));
      else
        diagonals[0]->matrix->vmult(aux, dst.block(row_num));
      dst.block(row_num).equ (diagonals[0]->prefix, aux);
    }
  else
    {
      aux = 0.;
      for (size_type i=0; i<diagonals.size(); ++i)
        {
          m = diagonals[i];
          // First, divide by the current
          // factor, such that we can
          // multiply by it later.
          aux /= m->prefix;
          if (m->transpose)
            m->matrix->Tvmult_add(aux, dst.block(row_num));
          else
            m->matrix->vmult_add(aux, dst.block(row_num));
          aux *= m->prefix;
        }
      dst.block(row_num) = aux;
    }
}



template <typename number, typename BLOCK_VECTOR>
void
BlockTrianglePrecondition<number,BLOCK_VECTOR>::vmult_add (
  BLOCK_VECTOR &dst,
  const BLOCK_VECTOR &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  BLOCK_VECTOR aux;
  aux.reinit(dst);
  vmult(aux, src);
  dst.add(aux);
}



template <typename number, typename BLOCK_VECTOR>
void
BlockTrianglePrecondition<number,BLOCK_VECTOR>::vmult (
  BLOCK_VECTOR &dst,
  const BLOCK_VECTOR &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  dst.equ(1., src);

  if (backward)
    {
      for (unsigned int i=n_block_rows(); i>0;)
        do_row(dst, --i);
    }
  else
    {
      for (unsigned int i=0; i<n_block_rows(); ++i)
        do_row(dst, i);
    }

}

template <typename number, typename BLOCK_VECTOR>
void
BlockTrianglePrecondition<number,BLOCK_VECTOR>::Tvmult (
  BLOCK_VECTOR &,
  const BLOCK_VECTOR &) const
{
  Assert (false, ExcNotImplemented());
}


template <typename number, typename BLOCK_VECTOR>
void
BlockTrianglePrecondition<number,BLOCK_VECTOR>::Tvmult_add (
  BLOCK_VECTOR &,
  const BLOCK_VECTOR &) const
{
  Assert (false, ExcNotImplemented());
}

template class BlockMatrixArray<float>;
template class BlockMatrixArray<double>;
template class BlockTrianglePrecondition<float>;
template class BlockTrianglePrecondition<double>;

#ifdef DEAL_II_WITH_TRILINOS
template class BlockMatrixArray<float, TrilinosWrappers::MPI::BlockVector>;
template class BlockMatrixArray<double, TrilinosWrappers::MPI::BlockVector>;
template class BlockTrianglePrecondition<float, TrilinosWrappers::MPI::BlockVector>;
template class BlockTrianglePrecondition<double, TrilinosWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_WITH_PETSC
template class BlockMatrixArray<float, PETScWrappers::MPI::BlockVector>;
template class BlockMatrixArray<double, PETScWrappers::MPI::BlockVector>;
template class BlockTrianglePrecondition<float, PETScWrappers::MPI::BlockVector>;
template class BlockTrianglePrecondition<double, PETScWrappers::MPI::BlockVector>;
#endif

DEAL_II_NAMESPACE_CLOSE
