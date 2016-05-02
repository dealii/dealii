// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2015 by the deal.II authors
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

#ifndef dealii__relaxation_block_templates_h
#define dealii__relaxation_block_templates_h

#include <deal.II/lac/relaxation_block.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

template <typename MatrixType, typename inverse_type>
inline
RelaxationBlock<MatrixType,inverse_type>::AdditionalData::AdditionalData
(const double relaxation,
 const bool   invert_diagonal,
 const bool   same_diagonal)
  :
  relaxation(relaxation),
  invert_diagonal(invert_diagonal),
  same_diagonal(same_diagonal),
  inversion(PreconditionBlockBase<inverse_type>::gauss_jordan),
  threshold(0.)
{}


template <typename MatrixType, typename inverse_type>
inline
std::size_t
RelaxationBlock<MatrixType,inverse_type>::AdditionalData::memory_consumption() const
{
  std::size_t result = sizeof(*this)
                       + block_list.memory_consumption()
                       - sizeof(block_list);
  for (unsigned int i=0; i<order.size(); ++i)
    result += MemoryConsumption::memory_consumption(order[i]);
  return result;
}


template <typename MatrixType, typename inverse_type>
inline
void
RelaxationBlock<MatrixType,inverse_type>::initialize (const MatrixType     &M,
                                                      const AdditionalData &parameters)
{
  Assert (parameters.invert_diagonal, ExcNotImplemented());

  clear();
//  Assert (M.m() == M.n(), ExcNotQuadratic());
  A = &M;
  additional_data = &parameters;
  this->inversion = parameters.inversion;

  this->reinit(additional_data->block_list.n_rows(), 0, additional_data->same_diagonal,
               additional_data->inversion);

  if (additional_data->invert_diagonal)
    invert_diagblocks();
}


template <typename MatrixType, typename inverse_type>
inline
void
RelaxationBlock<MatrixType,inverse_type>::clear ()
{
  A = 0;
  additional_data = 0;
  PreconditionBlockBase<inverse_type>::clear ();
}


template <typename MatrixType, typename inverse_type>
inline
void
RelaxationBlock<MatrixType,inverse_type>::invert_diagblocks ()
{
  const MatrixType &M=*A;
  FullMatrix<inverse_type> M_cell;

  if (this->same_diagonal())
    {
      Assert(false, ExcNotImplemented());
    }
  else
    {
      for (size_type block=0; block<additional_data->block_list.n_rows(); ++block)
        {
          const size_type bs = additional_data->block_list.row_length(block);
          M_cell.reinit(bs, bs);

          // Copy rows for this block
          // into the matrix for the
          // diagonal block
          SparsityPattern::iterator row
            = additional_data->block_list.begin(block);
          for (size_type row_cell=0; row_cell<bs; ++row_cell, ++row)
            {
//TODO:[GK] Optimize here
              for (typename MatrixType::const_iterator entry = M.begin(row->column());
                   entry != M.end(row->column()); ++entry)
                {
                  const size_type column = entry->column();
                  const size_type col_cell = additional_data->block_list.row_position(block, column);
                  if (col_cell != numbers::invalid_size_type)
                    M_cell(row_cell, col_cell) = entry->value();
                }
            }
          // Now M_cell contains the
          // diagonal block. Now
          // store it and its
          // inverse, if so requested.
          if (this->store_diagonals())
            {
              this->diagonal(block).reinit(bs, bs);
              this->diagonal(block) = M_cell;
            }
          switch (this->inversion)
            {
            case PreconditionBlockBase<inverse_type>::gauss_jordan:
              this->inverse(block).reinit(bs, bs);
              this->inverse(block).invert(M_cell);
              break;
            case PreconditionBlockBase<inverse_type>::householder:
              this->inverse_householder(block).initialize(M_cell);
              break;
            case PreconditionBlockBase<inverse_type>::svd:
              this->inverse_svd(block).reinit(bs, bs);
              this->inverse_svd(block) = M_cell;
              this->inverse_svd(block).compute_inverse_svd(additional_data->threshold);
              break;
            default:
              Assert(false, ExcNotImplemented());
            }
        }
    }
  this->inverses_computed(true);
}


template <typename MatrixType, typename inverse_type>
template <typename number2>
inline
void
RelaxationBlock<MatrixType,inverse_type>::do_step (Vector<number2>       &dst,
                                                   const Vector<number2> &prev,
                                                   const Vector<number2> &src,
                                                   const bool             backward) const
{
  Assert (additional_data->invert_diagonal, ExcNotImplemented());

  const MatrixType &M=*this->A;
  Vector<number2> b_cell, x_cell;

  const bool permutation_empty = additional_data->order.size() == 0;
  const unsigned int n_permutations = (permutation_empty)
                                      ? 1U : additional_data->order.size();
  const size_type n_blocks = additional_data->block_list.n_rows();

  if (!permutation_empty)
    for (unsigned int i=0; i<additional_data->order.size(); ++i)
      AssertDimension(additional_data->order[i].size(), this->size());

  for (unsigned int perm=0; perm<n_permutations; ++perm)
    {
      for (unsigned int bi=0; bi<n_blocks; ++bi)
        {
          const unsigned int raw_block = backward ? (n_blocks - bi - 1) : bi;
          const unsigned int block = permutation_empty
                                     ? raw_block
                                     : (backward
                                        ? (additional_data->order[n_permutations-1-perm][raw_block])
                                        : (additional_data->order[perm][raw_block]));

          const size_type bs = additional_data->block_list.row_length(block);

          b_cell.reinit(bs);
          x_cell.reinit(bs);
          // Collect off-diagonal parts
          SparsityPattern::iterator row = additional_data->block_list.begin(block);
          for (size_type row_cell=0; row_cell<bs; ++row_cell, ++row)
            {
              b_cell(row_cell) = src(row->column());
              for (typename MatrixType::const_iterator entry = M.begin(row->column());
                   entry != M.end(row->column()); ++entry)
                b_cell(row_cell) -= entry->value() * prev(entry->column());
            }
          // Apply inverse diagonal
          this->inverse_vmult(block, x_cell, b_cell);
#ifdef DEBUG
          for (unsigned int i=0; i<x_cell.size(); ++i)
            {
              AssertIsFinite(x_cell(i));
            }
#endif
          // Store in result vector
          row=additional_data->block_list.begin(block);
          for (size_type row_cell=0; row_cell<bs; ++row_cell, ++row)
            dst(row->column()) += additional_data->relaxation * x_cell(row_cell);
        }
    }
}


//----------------------------------------------------------------------//

template <typename MatrixType, typename inverse_type>
template <typename number2>
void RelaxationBlockJacobi<MatrixType,inverse_type>::step
(Vector<number2>       &dst,
 const Vector<number2> &src) const
{
  GrowingVectorMemory<Vector<number2> > mem;
  typename VectorMemory<Vector<number2> >::Pointer aux = mem;
  aux->reinit(dst, false);
  *aux = dst;
  this->do_step(dst, *aux, src, false);
}


template <typename MatrixType, typename inverse_type>
template <typename number2>
void RelaxationBlockJacobi<MatrixType,inverse_type>::Tstep
(Vector<number2>       &dst,
 const Vector<number2> &src) const
{
  GrowingVectorMemory<Vector<number2> > mem;
  typename VectorMemory<Vector<number2> >::Pointer aux = mem;
  aux->reinit(dst, false);
  *aux = dst;
  this->do_step(dst, *aux, src, true);
}


//----------------------------------------------------------------------//

template <typename MatrixType, typename inverse_type>
template <typename number2>
void RelaxationBlockSOR<MatrixType,inverse_type>::step
(Vector<number2> &dst,
 const Vector<number2> &src) const
{
  this->do_step(dst, dst, src, false);
}


template <typename MatrixType, typename inverse_type>
template <typename number2>
void RelaxationBlockSOR<MatrixType,inverse_type>::Tstep
(Vector<number2>       &dst,
 const Vector<number2> &src) const
{
  this->do_step(dst, dst, src, true);
}


//----------------------------------------------------------------------//

template <typename MatrixType, typename inverse_type>
template <typename number2>
void RelaxationBlockSSOR<MatrixType,inverse_type>::step
(Vector<number2>       &dst,
 const Vector<number2> &src) const
{
  this->do_step(dst, dst, src, false);
  this->do_step(dst, dst, src, true);
}


template <typename MatrixType, typename inverse_type>
template <typename number2>
void RelaxationBlockSSOR<MatrixType,inverse_type>::Tstep
(Vector<number2>       &dst,
 const Vector<number2> &src) const
{
  this->do_step(dst, dst, src, true);
  this->do_step(dst, dst, src, false);
}



DEAL_II_NAMESPACE_CLOSE


#endif
