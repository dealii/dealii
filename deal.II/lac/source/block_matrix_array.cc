//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/block_matrix_array.h>
#include <lac/vector.h>
#include <lac/block_vector.h>


template <typename number>
BlockMatrixArray<number>::Entry::Entry (const Entry& e)
		:
		row(e.row),
		col(e.col),
		prefix(e.prefix),
		transpose(e.transpose),
		matrix(e.matrix)
{
  Entry& e2 = const_cast<Entry&>(e);
  e2.matrix = 0;
}



template <typename number>
BlockMatrixArray<number>::Entry::~Entry ()
{
  if (matrix)
    delete matrix;
}



template <typename number>
BlockMatrixArray<number>::BlockMatrixArray (
  const unsigned int n_block_rows,
  const unsigned int n_block_cols,
  VectorMemory<Vector<number> >& mem)
		: block_rows (n_block_rows),
		  block_cols (n_block_cols),
		  mem(&mem, typeid(*this).name())
{}



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
BlockMatrixArray<number>::vmult_add (BlockVector<number>& dst,
				     const BlockVector<number>& src) const
{
  Assert (dst.n_blocks() == block_rows,
	  ExcDimensionMismatch(dst.n_blocks(), block_rows));
  Assert (src.n_blocks() == block_cols,
	  ExcDimensionMismatch(src.n_blocks(), block_cols));

  Vector<number>* p_aux = mem->alloc();
  Vector<number>& aux = *p_aux;
  
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
  mem->free(p_aux);
}




template <typename number>
void
BlockMatrixArray<number>::vmult (BlockVector<number>& dst,
				 const BlockVector<number>& src) const
{
  dst = 0.;
  vmult_add (dst, src);
}




template <typename number>
void
BlockMatrixArray<number>::Tvmult_add (BlockVector<number>& dst,
				      const BlockVector<number>& src) const
{
  Assert (dst.n_blocks() == block_cols,
	  ExcDimensionMismatch(dst.n_blocks(), block_cols));
  Assert (src.n_blocks() == block_rows,
	  ExcDimensionMismatch(src.n_blocks(), block_rows));

  typename std::vector<Entry>::const_iterator m = entries.begin();
  typename std::vector<Entry>::const_iterator end = entries.end();
  
  Vector<number>* p_aux = mem->alloc();
  Vector<number>& aux = *p_aux;
  
  for (; m != end ; ++m)
    {
      aux.reinit(dst.block(m->col));
      if (m->transpose)
	m->matrix->vmult(aux, src.block(m->row));
      else
	m->matrix->Tvmult(aux, src.block(m->row));
      dst.block(m->col).add (m->prefix, aux);
    }
  mem->free(p_aux);
}



template <typename number>
void
BlockMatrixArray<number>::Tvmult (BlockVector<number>& dst,
				  const BlockVector<number>& src) const
{
  dst = 0.;
  Tvmult_add (dst, src);
}




template <typename number>
number
BlockMatrixArray<number>::matrix_scalar_product (
  const BlockVector<number>& u,
  const BlockVector<number>& v) const
{
  Assert (u.n_blocks() == block_rows,
	  ExcDimensionMismatch(u.n_blocks(), block_rows));
  Assert (v.n_blocks() == block_cols,
	  ExcDimensionMismatch(v.n_blocks(), block_cols));

  Vector<number>* p_aux = mem->alloc();
  Vector<number>& aux = *p_aux;
  
  typename std::vector<Entry>::const_iterator m;
  typename std::vector<Entry>::const_iterator end = entries.end();

  number result = 0.;
  
  for (unsigned int i=0;i<block_rows;++i)
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
  mem->free(p_aux);

  return result;
}



template <typename number>
number
BlockMatrixArray<number>::matrix_norm_square (
  const BlockVector<number>& u) const
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
BlockTrianglePrecondition<number>::BlockTrianglePrecondition(
  unsigned int block_rows,
  VectorMemory<Vector<number> >& mem,
  bool backward)
		:
		BlockMatrixArray<number> (block_rows, block_rows, mem),
		backward(backward)
{}


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
  BlockVector<number>& dst,
  unsigned int row_num) const
{
  typename std::vector<typename BlockMatrixArray<number>::Entry>::const_iterator
    m = this->entries.begin();
  typename std::vector<typename BlockMatrixArray<number>::Entry>::const_iterator
    end = this->entries.end();
  typename std::vector<typename BlockMatrixArray<number>::Entry>::const_iterator
    diagonal = end;  
  
  Vector<number>* p_aux = this->mem->alloc();
  Vector<number>& aux = *p_aux;
  
  aux.reinit(dst.block(row_num), true);

  for (; m != end ; ++m)
    {
      const unsigned int i=m->row;
      if (i != row_num)
	continue;
      const unsigned int j=m->col;
      if (((j > i) && !backward) || ((j < i) && backward))
	continue;
      if (j == i)
	{
	  Assert (diagonal == end, ExcMultipleDiagonal(j));
	  diagonal = m;
	} else {
	  if (m->transpose)
	    m->matrix->Tvmult(aux, dst.block(j));
	  else
	    m->matrix->vmult(aux, dst.block(j));
	  dst.block(i).add (-1 * m->prefix, aux);
	}
    }
  if (diagonal->transpose)
    diagonal->matrix->Tvmult(aux, dst.block(row_num));
  else
    diagonal->matrix->vmult(aux, dst.block(row_num));
  dst.block(row_num).equ (diagonal->prefix, aux);
  
  this->mem->free(p_aux);
}



template <typename number>
void
BlockTrianglePrecondition<number>::vmult_add (
  BlockVector<number>& dst,
  const BlockVector<number>& src) const
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
  BlockVector<number>& dst,
  const BlockVector<number>& src) const
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
    } else {
      for (unsigned int i=0; i<n_block_rows(); ++i)
	do_row(dst, i);
    }
  
}

template <typename number>
void
BlockTrianglePrecondition<number>::Tvmult (
  BlockVector<number>&,
  const BlockVector<number>&) const
{
  Assert (false, ExcNotImplemented());
}


template <typename number>
void
BlockTrianglePrecondition<number>::Tvmult_add (
  BlockVector<number>&,
  const BlockVector<number>&) const
{
  Assert (false, ExcNotImplemented());
}

template class BlockMatrixArray<float>;
template class BlockMatrixArray<double>;
template class BlockTrianglePrecondition<float>;
template class BlockTrianglePrecondition<double>;
