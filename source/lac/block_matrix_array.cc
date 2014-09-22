// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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

DEAL_II_NAMESPACE_OPEN


template <typename number>
BlockMatrixArray<number>::Entry::Entry (const Entry &e)
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



template <typename number>
BlockMatrixArray<number>::Entry::~Entry ()
{
  if (matrix)
    delete matrix;
}



template <typename number>
BlockMatrixArray<number>::BlockMatrixArray ()
  : block_rows (0),
    block_cols (0)
{}



template <typename number>
BlockMatrixArray<number>::BlockMatrixArray (
  const unsigned int n_block_rows,
  const unsigned int n_block_cols)
  : block_rows (n_block_rows),
    block_cols (n_block_cols)
{}


template <typename number>
BlockMatrixArray<number>::BlockMatrixArray (
  const unsigned int n_block_rows,
  const unsigned int n_block_cols,
  VectorMemory<Vector<number> > &)
  : block_rows (n_block_rows),
    block_cols (n_block_cols)
{}


template <typename number>
void
BlockMatrixArray<number>::initialize (
  const unsigned int n_block_rows,
  const unsigned int n_block_cols,
  VectorMemory<Vector<number> > &)
{
  block_rows = n_block_rows;
  block_cols = n_block_cols;
}


template <typename number>
void
BlockMatrixArray<number>::initialize (
  const unsigned int n_block_rows,
  const unsigned int n_block_cols)
{
  block_rows = n_block_rows;
  block_cols = n_block_cols;
}



template <typename number>
void
BlockMatrixArray<number>::reinit (
  const unsigned int n_block_rows,
  const unsigned int n_block_cols)
{
  clear();
  block_rows = n_block_rows;
  block_cols = n_block_cols;
}



template <typename number>
void
BlockMatrixArray<number>::clear ()
{
  entries.clear();
}


template <typename number>
void
BlockMatrixArray<number>::vmult_add (BlockVector<number> &dst,
                                     const BlockVector<number> &src) const
{
  GrowingVectorMemory<Vector<number> > mem;
  Assert (dst.n_blocks() == block_rows,
          ExcDimensionMismatch(dst.n_blocks(), block_rows));
  Assert (src.n_blocks() == block_cols,
          ExcDimensionMismatch(src.n_blocks(), block_cols));

  typename VectorMemory<Vector<number> >::Pointer p_aux(mem);
  Vector<number> &aux = *p_aux;

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




template <typename number>
void
BlockMatrixArray<number>::vmult (BlockVector<number> &dst,
                                 const BlockVector<number> &src) const
{
  dst = 0.;
  vmult_add (dst, src);
}




template <typename number>
void
BlockMatrixArray<number>::Tvmult_add (BlockVector<number> &dst,
                                      const BlockVector<number> &src) const
{
  GrowingVectorMemory<Vector<number> > mem;
  Assert (dst.n_blocks() == block_cols,
          ExcDimensionMismatch(dst.n_blocks(), block_cols));
  Assert (src.n_blocks() == block_rows,
          ExcDimensionMismatch(src.n_blocks(), block_rows));

  typename std::vector<Entry>::const_iterator m = entries.begin();
  typename std::vector<Entry>::const_iterator end = entries.end();

  typename VectorMemory<Vector<number> >::Pointer p_aux(mem);
  Vector<number> &aux = *p_aux;

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



template <typename number>
void
BlockMatrixArray<number>::Tvmult (BlockVector<number> &dst,
                                  const BlockVector<number> &src) const
{
  dst = 0.;
  Tvmult_add (dst, src);
}




template <typename number>
number
BlockMatrixArray<number>::matrix_scalar_product (
  const BlockVector<number> &u,
  const BlockVector<number> &v) const
{
  GrowingVectorMemory<Vector<number> > mem;
  Assert (u.n_blocks() == block_rows,
          ExcDimensionMismatch(u.n_blocks(), block_rows));
  Assert (v.n_blocks() == block_cols,
          ExcDimensionMismatch(v.n_blocks(), block_cols));

  typename VectorMemory<Vector<number> >::Pointer p_aux(mem);
  Vector<number> &aux = *p_aux;

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



template <typename number>
number
BlockMatrixArray<number>::matrix_norm_square (
  const BlockVector<number> &u) const
{
  return matrix_scalar_product(u,u);
}



template <typename number>
unsigned int
BlockMatrixArray<number>::n_block_rows () const
{
  return block_rows;
}



template <typename number>
unsigned int
BlockMatrixArray<number>::n_block_cols () const
{
  return block_cols;
}



//---------------------------------------------------------------------------

template <typename number>
BlockTrianglePrecondition<number>::BlockTrianglePrecondition()
  : BlockMatrixArray<number> (),
    backward(false)
{}


template <typename number>
BlockTrianglePrecondition<number>::BlockTrianglePrecondition(
  const unsigned int block_rows,
  VectorMemory<Vector<number> > &,
  const bool backward)
  :
  BlockMatrixArray<number> (block_rows, block_rows),
  backward(backward)
{}


template <typename number>
BlockTrianglePrecondition<number>::BlockTrianglePrecondition(
  const unsigned int block_rows)
  :
  BlockMatrixArray<number> (block_rows, block_rows),
  backward(false)
{}


template <typename number>
void
BlockTrianglePrecondition<number>::initialize(
  const unsigned int n_block_rows,
  VectorMemory<Vector<number> > &,
  const bool backward)
{
  BlockMatrixArray<number>::initialize(n_block_rows, n_block_rows);
  this->backward = backward;
}


template <typename number>
void
BlockTrianglePrecondition<number>::reinit (
  const unsigned int n)
{
  BlockMatrixArray<number>::reinit(n,n);
}


template <typename number>
void
BlockTrianglePrecondition<number>::do_row (
  BlockVector<number> &dst,
  size_type row_num) const
{
  GrowingVectorMemory<Vector<number> > mem;
  typename std::vector<typename BlockMatrixArray<number>::Entry>::const_iterator
  m = this->entries.begin();
  typename std::vector<typename BlockMatrixArray<number>::Entry>::const_iterator
  end = this->entries.end();
  std::vector<typename std::vector<typename BlockMatrixArray<number>::Entry>::const_iterator>
  diagonals;

  typename VectorMemory<Vector<number> >::Pointer p_aux(mem);
  Vector<number> &aux = *p_aux;

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



template <typename number>
void
BlockTrianglePrecondition<number>::vmult_add (
  BlockVector<number> &dst,
  const BlockVector<number> &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  BlockVector<number> aux;
  aux.reinit(dst);
  vmult(aux, src);
  dst.add(aux);
}



template <typename number>
void
BlockTrianglePrecondition<number>::vmult (
  BlockVector<number> &dst,
  const BlockVector<number> &src) const
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

template <typename number>
void
BlockTrianglePrecondition<number>::Tvmult (
  BlockVector<number> &,
  const BlockVector<number> &) const
{
  Assert (false, ExcNotImplemented());
}


template <typename number>
void
BlockTrianglePrecondition<number>::Tvmult_add (
  BlockVector<number> &,
  const BlockVector<number> &) const
{
  Assert (false, ExcNotImplemented());
}

template class BlockMatrixArray<float>;
template class BlockMatrixArray<double>;
template class BlockTrianglePrecondition<float>;
template class BlockTrianglePrecondition<double>;

DEAL_II_NAMESPACE_CLOSE
