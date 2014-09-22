// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__precondition_block_templates_h
#define __deal2__precondition_block_templates_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/householder.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

DEAL_II_NAMESPACE_OPEN


template<class MATRIX, typename inverse_type>
PreconditionBlock<MATRIX, inverse_type>::AdditionalData::
AdditionalData (const size_type block_size,
                const double relaxation,
                const bool invert_diagonal,
                const bool same_diagonal)
  :
  relaxation (relaxation),
  block_size(block_size),
  invert_diagonal(invert_diagonal),
  same_diagonal(same_diagonal),
  inversion(PreconditionBlockBase<inverse_type>::gauss_jordan),
  threshold(0.)
{}


template <typename number>
PreconditionBlockBase<number>::~PreconditionBlockBase ()
{}


template <class MATRIX, typename inverse_type>
PreconditionBlock<MATRIX,inverse_type>::PreconditionBlock (bool store)
  : PreconditionBlockBase<inverse_type>(store),
    blocksize(0),
    A(0, typeid(*this).name())
{}


template <class MATRIX, typename inverse_type>
PreconditionBlock<MATRIX,inverse_type>::~PreconditionBlock ()
{}


template <class MATRIX, typename inverse_type>
void PreconditionBlock<MATRIX,inverse_type>::clear ()
{
  PreconditionBlockBase<inverse_type>::clear();
  blocksize     = 0;
  A = 0;
}


template <class MATRIX, typename inverse_type>
void PreconditionBlock<MATRIX,inverse_type>::initialize (
  const MATRIX &M,
  const AdditionalData parameters)
{
  const size_type bsize = parameters.block_size;

  clear();
  Assert (M.m() == M.n(), ExcNotQuadratic());
  A = &M;
  Assert (bsize>0, ExcIndexRange(bsize, 1, M.m()));
  Assert (A->m()%bsize==0, ExcWrongBlockSize(bsize, A->m()));
  blocksize=bsize;
  relaxation = parameters.relaxation;
  const unsigned int nblocks = A->m()/bsize;
  this->reinit(nblocks, blocksize, parameters.same_diagonal,
               parameters.inversion);

  if (parameters.invert_diagonal)
    {
      if (permutation.size() == M.m())
        invert_permuted_diagblocks(permutation, inverse_permutation);
      else
        invert_diagblocks();
    }
}


template <class MATRIX, typename inverse_type>
void PreconditionBlock<MATRIX,inverse_type>::initialize (
  const MATRIX &M,
  const std::vector<size_type> &permutation,
  const std::vector<size_type> &inverse_permutation,
  const AdditionalData parameters)
{
  set_permutation(permutation, inverse_permutation);
  initialize(M, parameters);
}

template <class MATRIX, typename inverse_type>
void PreconditionBlock<MATRIX,inverse_type>::invert_permuted_diagblocks(
  const std::vector<size_type> &permutation,
  const std::vector<size_type> &inverse_permutation)
{
  Assert (A!=0, ExcNotInitialized());
  Assert (blocksize!=0, ExcNotInitialized());

  const MATRIX &M=*A;
  Assert (this->inverses_ready()==0, ExcInverseMatricesAlreadyExist());
  AssertDimension (permutation.size(), M.m());
  AssertDimension (inverse_permutation.size(), M.m());

  FullMatrix<inverse_type> M_cell(blocksize);

  if (this->same_diagonal())
    {
      deallog << "PreconditionBlock uses only one diagonal block" << std::endl;

      for (size_type row_cell=0; row_cell<blocksize; ++row_cell)
        {
          typename MATRIX::const_iterator entry = M.begin(row_cell);
          const typename MATRIX::const_iterator row_end = M.end(row_cell);
          while (entry != row_end)
            {
              if (entry->column() < blocksize)
                M_cell(row_cell, entry->column()) = entry->value();
              ++entry;
            }
        }
      if (this->store_diagonals())
        this->diagonal(0) = M_cell;
      switch (this->inversion)
        {
        case PreconditionBlockBase<inverse_type>::gauss_jordan:
          this->inverse(0).invert(M_cell);
          break;
        case PreconditionBlockBase<inverse_type>::householder:
          this->inverse_householder(0).initialize(M_cell);
          break;
        case PreconditionBlockBase<inverse_type>::svd:
          this->inverse_svd(0) = M_cell;
          this->inverse_svd(0).compute_inverse_svd(0.);
          break;
        default:
          Assert(false, ExcNotImplemented());

        }
    }
  else
    {
      // cell_row, cell_column are the
      // numbering of the blocks (cells).
      // row_cell, column_cell are the local
      // numbering of the unknowns in the
      // blocks.
      // row, column are the global numbering
      // of the unkowns.
      M_cell = 0;

      for (unsigned int cell=0; cell<this->size(); ++cell)
        {
          const size_type cell_start = cell*blocksize;
          for (size_type row_cell=0; row_cell<blocksize; ++row_cell)
            {
              const size_type urow = row_cell + cell_start;

              const size_type row = permutation[urow];

              typename MATRIX::const_iterator entry = M.begin(row);
              const typename MATRIX::const_iterator row_end = M.end(row);

              for (; entry != row_end; ++entry)
                {
                  //if (entry->column()<cell_start)
                  if (inverse_permutation[entry->column()]<cell_start)
                    continue;

                  const size_type column_cell = inverse_permutation[entry->column()]-cell_start;
                  if (column_cell >= blocksize)
                    continue;
                  M_cell(row_cell, column_cell) = entry->value();
                }
            }

          if (this->store_diagonals())
            this->diagonal(cell) = M_cell;
          switch (this->inversion)
            {
            case PreconditionBlockBase<inverse_type>::gauss_jordan:
              this->inverse(cell).invert(M_cell);
              break;
            case PreconditionBlockBase<inverse_type>::householder:
              this->inverse_householder(cell).initialize(M_cell);
              break;
            case PreconditionBlockBase<inverse_type>::svd:
              this->inverse_svd(cell) = M_cell;
              this->inverse_svd(cell).compute_inverse_svd(0.);
              break;
            default:
              Assert(false, ExcNotImplemented());

            }
        }
    }
  this->inverses_computed(true);
}



template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlock<MATRIX,inverse_type>::forward_step (
  Vector<number2>       &dst,
  const Vector<number2> &prev,
  const Vector<number2> &src,
  const bool transpose_diagonal) const
{
  Assert (this->A!=0, ExcNotInitialized());

  const MATRIX &M=*this->A;

  if (permutation.size() != 0)
    Assert (permutation.size() == M.m() || permutation.size() == this->size(),
            ExcMessage("Permutation vector size must be equal to either the number of blocks or the dimension of the system"));

  const bool permuted = (permutation.size() == M.m());
  const bool cell_permuted = (permutation.size() == this->size());

//   deallog << "Permutation " << permutation.size();
//   if (permuted) deallog << " point";
//   if (cell_permuted) deallog << " block";
//   deallog << std::endl;

  Vector<number2> b_cell(this->blocksize), x_cell(this->blocksize);

  // cell_row, cell_column are the
  // numbering of the blocks (cells).
  // row_cell, column_cell are the local
  // numbering of the unknowns in the
  // blocks.
  // row, column are the global numbering
  // of the unkowns.
  size_type row, row_cell;
  number2 b_cell_row;
  // The diagonal block if the
  // inverses were not precomputed
  FullMatrix<number> M_cell(this->blocksize);

  // Loop over all blocks
  for (unsigned int rawcell=0; rawcell < this->size(); ++rawcell)
    {
      const unsigned int cell = cell_permuted ? permutation[rawcell] : rawcell;
      const size_type block_start = cell*this->blocksize;
      const size_type permuted_block_start = permuted
                                             ? permutation[block_start]
                                             : block_start;

//       deallog << std::endl << cell << '-' << block_start
//            << '-' << permuted_block_start << (permuted ? 't' : 'f') << '\t';

      for (row = permuted_block_start, row_cell = 0;
           row_cell < this->blocksize;
           ++row_cell, ++row)
        {
//        deallog << ' ' << row;
          const typename MATRIX::const_iterator row_end = M.end(row);
          typename MATRIX::const_iterator entry = M.begin(row);

          b_cell_row=src(row);
          for (; entry != row_end; ++entry)
            {
              const size_type column = entry->column();
              const size_type inverse_permuted_column = permuted
                                                        ? inverse_permutation[column]
                                                        : column;
              b_cell_row -= entry->value() * prev(column);
//TODO:[GK] Find out if this is really once column and once permuted
              if (!this->inverses_ready()
                  && inverse_permuted_column >= block_start
                  && inverse_permuted_column < block_start + this->blocksize)
                {
                  const size_type column_cell = column - block_start;
                  if (transpose_diagonal)
                    M_cell(column_cell, row_cell) = entry->value();
                  else
                    M_cell(row_cell, column_cell) = entry->value();
                }
            }
          b_cell(row_cell)=b_cell_row;
        }
      if (this->inverses_ready())
        {
          if (transpose_diagonal)
            this->inverse_Tvmult(cell, x_cell, b_cell);
          else
            this->inverse_vmult(cell, x_cell, b_cell);
        }
      else
        {
          Householder<number> house(M_cell);
          house.least_squares(x_cell,b_cell);
        }

      // distribute x_cell to dst
      for (row=permuted_block_start, row_cell=0;
           row_cell<this->blocksize;
           ++row_cell, ++row)
        dst(row) = prev(row) + this->relaxation*x_cell(row_cell);
    }
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlock<MATRIX,inverse_type>::backward_step (
  Vector<number2>       &dst,
  const Vector<number2> &prev,
  const Vector<number2> &src,
  const bool transpose_diagonal) const
{
  Assert (this->A!=0, ExcNotInitialized());

  const MATRIX &M=*this->A;

  if (permutation.size() != 0)
    Assert (permutation.size() == M.m() || permutation.size() == this->size(),
            ExcMessage("Permutation vector size must be equal to either the number of blocks or the dimension of the system"));

  const bool permuted = (permutation.size() == M.m());
  const bool cell_permuted = (permutation.size() == this->size());

  Vector<number2> b_cell(this->blocksize), x_cell(this->blocksize);

  // cell_row, cell_column are the
  // numbering of the blocks (cells).
  // row_cell, column_cell are the local
  // numbering of the unknowns in the
  // blocks.
  // row, column are the global numbering
  // of the unkowns.
  size_type row, row_cell;
  number2 b_cell_row;

  FullMatrix<number> M_cell(this->blocksize);
  for (unsigned int rawcell=this->size(); rawcell!=0 ;)
    {
      --rawcell;
      const unsigned int cell = cell_permuted ? permutation[rawcell] : rawcell;
      const size_type block_start = cell*this->blocksize;
      const size_type block_end = block_start + this->blocksize;
      const size_type permuted_block_start = permuted
                                             ? permutation[block_start]
                                             : block_start;
      for (row = permuted_block_start, row_cell = 0;
           row_cell<this->blocksize;
           ++row_cell, ++row)
        {
          const typename MATRIX::const_iterator row_end = M.end(row);
          typename MATRIX::const_iterator entry = M.begin(row);

          b_cell_row=src(row);
          for (; entry != row_end; ++entry)
            {
              const size_type column = entry->column();
              const size_type inverse_permuted_column = permuted
                                                        ? inverse_permutation[column]
                                                        : column;
              b_cell_row -= entry->value() * prev(column);
              if (!this->inverses_ready()
                  && inverse_permuted_column < block_end
                  && column >= block_start)
                {
                  const size_type column_cell = column - block_start;
                  // We need the
                  // transpose of the
                  // diagonal block,
                  // so we switch row
                  // and column
                  // indices
                  if (transpose_diagonal)
                    M_cell(column_cell, row_cell) = entry->value();
                  else
                    M_cell(row_cell, column_cell) = entry->value();
                }
            }
          b_cell(row_cell)=b_cell_row;
        }
      if (this->inverses_ready())
        {
          if (transpose_diagonal)
            this->inverse_Tvmult(cell, x_cell, b_cell);
          else
            this->inverse_vmult(cell, x_cell, b_cell);
        }
      else
        {
          Householder<number> house(M_cell);
          house.least_squares(x_cell,b_cell);
        }


      // distribute x_cell to dst
      for (row=permuted_block_start, row_cell=0;
           row_cell<this->blocksize;
           ++row_cell, ++row)
        dst(row) = prev(row) + this->relaxation*x_cell(row_cell);
    }
}


template <class MATRIX, typename inverse_type>
typename PreconditionBlock<MATRIX,inverse_type>::size_type
PreconditionBlock<MATRIX,inverse_type>::block_size() const
{
  return blocksize;
}


template <class MATRIX, typename inverse_type>
void PreconditionBlock<MATRIX,inverse_type>::invert_diagblocks()
{
  Assert (A!=0, ExcNotInitialized());
  Assert (blocksize!=0, ExcNotInitialized());

  const MATRIX &M=*A;
  Assert (this->inverses_ready()==0, ExcInverseMatricesAlreadyExist());

  FullMatrix<inverse_type> M_cell(blocksize);

  if (this->same_diagonal())
    {
      deallog << "PreconditionBlock uses only one diagonal block" << std::endl;
      for (size_type row_cell=0; row_cell<blocksize; ++row_cell)
        {
          typename MATRIX::const_iterator entry = M.begin(row_cell);
          const typename MATRIX::const_iterator row_end = M.end(row_cell);
          while (entry != row_end)
            {
              if (entry->column() < blocksize)
                M_cell(row_cell, entry->column()) = entry->value();
              ++entry;
            }
        }
      if (this->store_diagonals())
        this->diagonal(0) = M_cell;
      switch (this->inversion)
        {
        case PreconditionBlockBase<inverse_type>::gauss_jordan:
          this->inverse(0).invert(M_cell);
          break;
        case PreconditionBlockBase<inverse_type>::householder:
          this->inverse_householder(0).initialize(M_cell);
          break;
        case PreconditionBlockBase<inverse_type>::svd:
          this->inverse_svd(0) = M_cell;
          this->inverse_svd(0).compute_inverse_svd(1.e-12);
          break;
        default:
          Assert(false, ExcNotImplemented());
        }
    }
  else
    {
      M_cell = 0;

      for (unsigned int cell=0; cell<this->size(); ++cell)
        {
          const size_type cell_start = cell*blocksize;
          for (size_type row_cell=0; row_cell<blocksize; ++row_cell)
            {
              const size_type row = row_cell + cell_start;
              typename MATRIX::const_iterator entry = M.begin(row);
              const typename MATRIX::const_iterator row_end = M.end(row);

              for (; entry != row_end; ++entry)
                {
                  if (entry->column()<cell_start)
                    continue;

                  const size_type column_cell = entry->column()-cell_start;
                  if (column_cell >= blocksize)
                    continue;
                  M_cell(row_cell, column_cell) = entry->value();
                }
            }

          if (this->store_diagonals())
            this->diagonal(cell) = M_cell;
          switch (this->inversion)
            {
            case PreconditionBlockBase<inverse_type>::gauss_jordan:
              this->inverse(cell).invert(M_cell);
              break;
            case PreconditionBlockBase<inverse_type>::householder:
              this->inverse_householder(cell).initialize(M_cell);
              break;
            case PreconditionBlockBase<inverse_type>::svd:
              this->inverse_svd(cell) = M_cell;
              this->inverse_svd(cell).compute_inverse_svd(1.e-12);
              break;
            default:
              Assert(false, ExcNotImplemented());
            }
        }
    }
  this->inverses_computed(true);
}



template <class MATRIX, typename inverse_type>
void PreconditionBlock<MATRIX,inverse_type>::set_permutation (
  const std::vector<size_type> &p,
  const std::vector<size_type> &i)
{
  Assert (p.size() == i.size(), ExcDimensionMismatch(p.size(), i.size()));

  if (this->inverses_ready())
    {
      AssertDimension(p.size(), this->size());
    }

  permutation.resize(p.size());
  inverse_permutation.resize(p.size());
  for (unsigned int k=0; k<p.size(); ++k)
    {
      permutation[k] = p[k];
      inverse_permutation[k] = i[k];
    }
}


template <class MATRIX, typename inverse_type>
std::size_t
PreconditionBlock<MATRIX,inverse_type>::memory_consumption () const
{
  return (sizeof(*this)
          - sizeof(PreconditionBlockBase<inverse_type>)
          + PreconditionBlockBase<inverse_type>::memory_consumption());
}




/*--------------------- PreconditionBlockJacobi -----------------------*/


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockJacobi<MATRIX,inverse_type>
::do_vmult (Vector<number2>       &dst,
            const Vector<number2> &src,
            bool adding) const
{
  Assert(this->A!=0, ExcNotInitialized());

  const MATRIX &M=*this->A;

  Vector<number2> b_cell(this->blocksize), x_cell(this->blocksize);

  // cell_row, cell_column are the
  // numbering of the blocks (cells).
  // row_cell, column_cell are the local
  // numbering of the unknowns in the
  // blocks.
  // row, column are the global numbering
  // of the unkowns.
  size_type row, row_cell, begin_diag_block=0;

  if (!this->inverses_ready())
    {
      FullMatrix<number> M_cell(this->blocksize);
      for (unsigned int cell=0; cell < this->size(); ++cell)
        {
          for (row=cell*this->blocksize, row_cell=0;
               row_cell<this->blocksize;
               ++row_cell, ++row)
            {
              b_cell(row_cell)=src(row);
              for (size_type column_cell=0, column=cell*this->blocksize;
                   column_cell<this->blocksize; ++column_cell, ++column)
                M_cell(row_cell,column_cell)=M(row,column);
            }
          Householder<number> house(M_cell);
          house.least_squares(x_cell,b_cell);
          // distribute x_cell to dst
          for (row=cell*this->blocksize, row_cell=0;
               row_cell<this->blocksize;
               ++row_cell, ++row)
            if (adding)
              dst(row)+=x_cell(row_cell);
            else
              dst(row)=x_cell(row_cell);

          begin_diag_block+=this->blocksize;
        }
    }
  else
    for (unsigned int cell=0; cell < this->size(); ++cell)
      {
        for (row=cell*this->blocksize, row_cell=0;
             row_cell<this->blocksize;
             ++row_cell, ++row)
          {
            b_cell(row_cell)=src(row);
          }
        this->inverse_vmult(cell, x_cell, b_cell);
        // distribute x_cell to dst
        for (row=cell*this->blocksize, row_cell=0;
             row_cell<this->blocksize;
             ++row_cell, ++row)
          if (adding)
            dst(row)+=x_cell(row_cell);
          else
            dst(row)=x_cell(row_cell);

        begin_diag_block+=this->blocksize;
      }
  dst *= this->relaxation;
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockJacobi<MATRIX,inverse_type>
::vmult (Vector<number2>       &dst,
         const Vector<number2> &src) const
{
  do_vmult(dst, src, false);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockJacobi<MATRIX,inverse_type>
::Tvmult (Vector<number2>       &dst,
          const Vector<number2> &src) const
{
  do_vmult(dst, src, false);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockJacobi<MATRIX,inverse_type>
::vmult_add (Vector<number2>       &dst,
             const Vector<number2> &src) const
{
  do_vmult(dst, src, true);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockJacobi<MATRIX,inverse_type>
::Tvmult_add (Vector<number2>       &dst,
              const Vector<number2> &src) const
{
  do_vmult(dst, src, true);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockJacobi<MATRIX,inverse_type>
::step (Vector<number2>       &dst,
        const Vector<number2> &src) const
{
  GrowingVectorMemory<Vector<number2> > mem;
  typename VectorMemory<Vector<number2> >::Pointer aux(mem);
  aux->reinit(dst);

  this->forward_step(*aux, dst, src, false);
  dst = *aux;
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockJacobi<MATRIX,inverse_type>
::Tstep (Vector<number2>       &dst,
         const Vector<number2> &src) const
{
  GrowingVectorMemory<Vector<number2> > mem;
  typename VectorMemory<Vector<number2> >::Pointer aux(mem);
  aux->reinit(dst);

  this->backward_step(*aux, dst, src, true);
  dst = *aux;
}




/*--------------------- PreconditionBlockSOR -----------------------*/


template <class MATRIX, typename inverse_type>
PreconditionBlockSOR<MATRIX,inverse_type>::PreconditionBlockSOR ()
  : PreconditionBlock<MATRIX,inverse_type> (false)

{}

template <class MATRIX, typename inverse_type>
PreconditionBlockSOR<MATRIX,inverse_type>::PreconditionBlockSOR (bool store)
  : PreconditionBlock<MATRIX,inverse_type> (store)

{}

template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSOR<MATRIX,inverse_type>::forward (
  Vector<number2>       &dst,
  const Vector<number2> &src,
  const bool transpose_diagonal,
  const bool) const
{
  Assert (this->A!=0, ExcNotInitialized());

  const MATRIX &M=*this->A;
  const bool permuted = (this->permutation.size() != 0);
  if (permuted)
    {
      Assert (this->permutation.size() == M.m(), ExcDimensionMismatch(this->permutation.size(), M.m()));
    }

  Vector<number2> b_cell(this->blocksize), x_cell(this->blocksize);

  // cell_row, cell_column are the
  // numbering of the blocks (cells).
  // row_cell, column_cell are the local
  // numbering of the unknowns in the
  // blocks.
  // row, column are the global numbering
  // of the unkowns.
  size_type row, row_cell, block_start=0;
  number2 b_cell_row;
  // The diagonal block if the
  // inverses were not precomputed
  FullMatrix<number> M_cell(this->blocksize);

  for (unsigned int cell=0; cell < this->size(); ++cell)
    {
      const size_type permuted_block_start = permuted
                                             ? this->permutation[block_start]
                                             :block_start;

      for (row = permuted_block_start, row_cell = 0;
           row_cell < this->blocksize;
           ++row_cell, ++row)
        {
          const typename MATRIX::const_iterator row_end = M.end(row);
          typename MATRIX::const_iterator entry = M.begin(row);

          b_cell_row=src(row);
          for (; entry != row_end; ++entry)
            {
              const size_type column = entry->column();
              const size_type inverse_permuted_column = permuted
                                                        ? this->inverse_permutation[column]
                                                        : column;

              if (inverse_permuted_column < block_start)
                b_cell_row -= entry->value() * dst(column);
              else if (!this->inverses_ready() && column < block_start + this->blocksize)
                {
                  const size_type column_cell = column - block_start;
                  if (transpose_diagonal)
                    M_cell(column_cell, row_cell) = entry->value();
                  else
                    M_cell(row_cell, column_cell) = entry->value();
                }
            }
          b_cell(row_cell)=b_cell_row;
        }
      if (this->inverses_ready())
        {
          if (transpose_diagonal)
            this->inverse_Tvmult(cell, x_cell, b_cell);
          else
            this->inverse_vmult(cell, x_cell, b_cell);
        }
      else
        {
          Householder<number> house(M_cell);
          house.least_squares(x_cell,b_cell);
        }

      // distribute x_cell to dst
      for (row=permuted_block_start, row_cell=0;
           row_cell<this->blocksize;
           ++row_cell, ++row)
        dst(row)=this->relaxation*x_cell(row_cell);

      block_start+=this->blocksize;
    }
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSOR<MATRIX,inverse_type>::backward (
  Vector<number2>       &dst,
  const Vector<number2> &src,
  const bool transpose_diagonal,
  const bool) const
{
  Assert (this->A!=0, ExcNotInitialized());

  const MATRIX &M=*this->A;
  const bool permuted = (this->permutation.size() != 0);
  if (permuted)
    {
      Assert (this->permutation.size() == M.m(), ExcDimensionMismatch(this->permutation.size(), M.m()));
    }

  Vector<number2> b_cell(this->blocksize), x_cell(this->blocksize);

  // cell_row, cell_column are the
  // numbering of the blocks (cells).
  // row_cell, column_cell are the local
  // numbering of the unknowns in the
  // blocks.
  // row, column are the global numbering
  // of the unkowns.
  size_type row, row_cell;
  size_type block_end=this->blocksize * this->size();
  number2 b_cell_row;

  FullMatrix<number> M_cell(this->blocksize);
  for (unsigned int cell=this->size(); cell!=0 ;)
    {
      --cell;
      const size_type block_start = block_end - this->blocksize;
      // Collect upper triangle
      const size_type permuted_block_start = (this->permutation.size() != 0)
                                             ? this->permutation[block_start]
                                             : block_start;
      for (row = permuted_block_start, row_cell = 0;
           row_cell<this->blocksize;
           ++row_cell, ++row)
        {
          const typename MATRIX::const_iterator row_end = M.end(row);
          typename MATRIX::const_iterator entry = M.begin(row);

          b_cell_row=src(row);
          for (; entry != row_end; ++entry)
            {
              const size_type column = entry->column();
              const size_type inverse_permuted_column = permuted
                                                        ? this->inverse_permutation[column]
                                                        : column;
              if (inverse_permuted_column >= block_end)
                b_cell_row -= entry->value() * dst(column);
              else if (!this->inverses_ready() && column >= block_start)
                {
                  const size_type column_cell = column - block_start;
                  // We need the
                  // transpose of the
                  // diagonal block,
                  // so we switch row
                  // and column
                  // indices
                  if (transpose_diagonal)
                    M_cell(column_cell, row_cell) = entry->value();
                  else
                    M_cell(row_cell, column_cell) = entry->value();
                }
            }
          b_cell(row_cell)=b_cell_row;
        }
      if (this->inverses_ready())
        {
          if (transpose_diagonal)
            this->inverse_Tvmult(cell, x_cell, b_cell);
          else
            this->inverse_vmult(cell, x_cell, b_cell);
        }
      else
        {
          Householder<number> house(M_cell);
          house.least_squares(x_cell,b_cell);
        }


      // distribute x_cell to dst
      for (row=permuted_block_start, row_cell=0;
           row_cell<this->blocksize;
           ++row_cell, ++row)
        dst(row)=this->relaxation*x_cell(row_cell);
      block_end = block_start;

    }
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSOR<MATRIX,inverse_type>
::vmult (Vector<number2>       &dst,
         const Vector<number2> &src) const
{
  forward(dst, src, false, false);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSOR<MATRIX,inverse_type>
::vmult_add (Vector<number2>       &dst,
             const Vector<number2> &src) const
{
  forward(dst, src, false, true);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSOR<MATRIX,inverse_type>
::Tvmult (Vector<number2>       &dst,
          const Vector<number2> &src) const
{
  backward(dst, src, true, false);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSOR<MATRIX,inverse_type>
::Tvmult_add (Vector<number2>       &dst,
              const Vector<number2> &src) const
{
  backward(dst, src, true, true);
}



template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSOR<MATRIX,inverse_type>
::step (Vector<number2>       &dst,
        const Vector<number2> &src) const
{
  this->forward_step(dst, dst, src, false);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSOR<MATRIX,inverse_type>
::Tstep (Vector<number2>       &dst,
         const Vector<number2> &src) const
{
  this->backward_step(dst, dst, src, true);
}




//---------------------------------------------------------------------------


template <class MATRIX, typename inverse_type>
PreconditionBlockSSOR<MATRIX,inverse_type>::PreconditionBlockSSOR ()
  : PreconditionBlockSOR<MATRIX,inverse_type> (true)

{}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSSOR<MATRIX,inverse_type>::vmult (Vector<number2>       &dst,
                                                        const Vector<number2> &src) const
{
  Vector<number2> help;
  help.reinit(dst);

  this->forward(help, src, false, false);

  Vector<inverse_type> cell_src(this->blocksize);
  Vector<inverse_type> cell_dst(this->blocksize);
  const double scaling = (2.-this->relaxation)/this->relaxation;

  // Multiply with diagonal blocks
  for (unsigned int cell=0; cell < this->size(); ++cell)
    {
      size_type row = cell*this->blocksize;

      for (size_type row_cell=0; row_cell<this->blocksize; ++row_cell)
        cell_src(row_cell)=help(row+row_cell);

      this->diagonal(cell).vmult(cell_dst, cell_src);

      for (size_type row_cell=0; row_cell<this->blocksize; ++row_cell)
        help(row+row_cell) = scaling * cell_dst(row_cell);
    }

  this->backward(dst, help, false, false);
}

template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSSOR<MATRIX,inverse_type>::Tvmult (Vector<number2>       &dst,
                                                         const Vector<number2> &src) const
{
  Vector<number2> help;
  help.reinit(dst);

  this->backward(help, src, true, false);

  Vector<inverse_type> cell_src(this->blocksize);
  Vector<inverse_type> cell_dst(this->blocksize);
  const double scaling = (2.-this->relaxation)/this->relaxation;

  // Multiply with diagonal blocks
  for (unsigned int cell=0; cell < this->size(); ++cell)
    {
      size_type row = cell*this->blocksize;

      for (size_type row_cell=0; row_cell<this->blocksize; ++row_cell)
        cell_src(row_cell)=help(row+row_cell);

      this->diagonal(cell).Tvmult(cell_dst, cell_src);

      for (size_type row_cell=0; row_cell<this->blocksize; ++row_cell)
        help(row+row_cell) = scaling * cell_dst(row_cell);
    }

  this->forward(dst, help, true, false);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSSOR<MATRIX,inverse_type>
::step (Vector<number2>       &dst,
        const Vector<number2> &src) const
{
  this->forward_step(dst, dst, src, false);
  this->backward_step(dst, dst, src, false);
}


template <class MATRIX, typename inverse_type>
template <typename number2>
void PreconditionBlockSSOR<MATRIX,inverse_type>
::Tstep (Vector<number2>       &dst,
         const Vector<number2> &src) const
{
  this->backward_step(dst, dst, src, true);
  this->forward_step(dst, dst, src, true);
}



DEAL_II_NAMESPACE_CLOSE

#endif
