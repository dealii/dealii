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
#include <deal.II/lac/trilinos_vector.h>

DEAL_II_NAMESPACE_OPEN

template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::AdditionalData::AdditionalData
(const double relaxation,
 const bool   invert_diagonal,
 const bool   same_diagonal,
 const typename PreconditionBlockBase<InverseNumberType>::Inversion inversion,
 const double threshold,
 VectorType *temp_ghost_vector)
  :
  relaxation(relaxation),
  invert_diagonal(invert_diagonal),
  same_diagonal(same_diagonal),
  inversion(inversion),
  threshold(threshold),
  temp_ghost_vector (temp_ghost_vector)
{}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline
std::size_t
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::AdditionalData::memory_consumption() const
{
  std::size_t result = sizeof(*this)
                       + block_list.memory_consumption()
                       - sizeof(block_list);
  for (unsigned int i=0; i<order.size(); ++i)
    result += MemoryConsumption::memory_consumption(order[i]);
  return result;
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline
void
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::initialize (const MatrixType     &M,
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


template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline
void
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::clear ()
{
  A = 0;
  additional_data = 0;
  PreconditionBlockBase<InverseNumberType>::clear ();
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline
void
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::invert_diagblocks ()
{
  const MatrixType &M=*A;
  FullMatrix<InverseNumberType> M_cell;

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
            case PreconditionBlockBase<InverseNumberType>::gauss_jordan:
              this->inverse(block).reinit(bs, bs);
              this->inverse(block).invert(M_cell);
              break;
            case PreconditionBlockBase<InverseNumberType>::householder:
              this->inverse_householder(block).initialize(M_cell);
              break;
            case PreconditionBlockBase<InverseNumberType>::svd:
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

namespace internal
{
  /**
   * Default implementation for serial vectors. Here we don't need to make a
   * copy into a ghosted vector, so just return a reference to @p prev.
   */
  template <class VectorType>
  const VectorType &
  prepare_ghost_vector(
    const VectorType &prev,
    VectorType *other)
  {
    // If the following Assertion triggers, you either set temp_ghost_vector
    // for a serial computation (don't!), or nobody implemented, instantiated, and
    // tested the parallel version for your vector type.
    Assert(other==NULL, ExcNotImplemented());
    (void)other;
    return prev;
  }

  /**
   * Specialization for Trilinos. Use the ghosted vector.
   */
  template <>
  const TrilinosWrappers::MPI::Vector &
  prepare_ghost_vector(
    const TrilinosWrappers::MPI::Vector &prev,
    TrilinosWrappers::MPI::Vector *other)
  {
    Assert(other!=NULL,
           ExcMessage("You need to provide a ghosted vector in RelaxationBlock::AdditionalData::temp_trilinos_ghost_vector."));
    Assert(other->size()==prev.size(), ExcInternalError());

    // import ghost values:
    *other = prev;
    return *other;
  }
} // end namespace internal

template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline
void
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::do_step (VectorType       &dst,
    const VectorType &prev,
    const VectorType &src,
    const bool             backward) const
{
  Assert (additional_data->invert_diagonal, ExcNotImplemented());

  const VectorType &ghosted_prev = internal::prepare_ghost_vector(prev, additional_data->temp_ghost_vector);

  const MatrixType &M=*this->A;
  Vector<typename VectorType::value_type> b_cell, x_cell;

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
                b_cell(row_cell) -= entry->value() * ghosted_prev(entry->column());
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
          row = additional_data->block_list.begin(block);
          for (size_type row_cell=0; row_cell<bs; ++row_cell, ++row)
            dst(row->column()) += additional_data->relaxation * x_cell(row_cell);
        }
    }
  dst.compress(dealii::VectorOperation::add);
}


//----------------------------------------------------------------------//

template <typename MatrixType, typename InverseNumberType, typename VectorType>
void RelaxationBlockJacobi<MatrixType, InverseNumberType, VectorType>::step
(VectorType       &dst,
 const VectorType &src) const
{
  GrowingVectorMemory<VectorType> mem;
  typename VectorMemory<VectorType>::Pointer aux = mem;
  aux->reinit(dst, false);
  *aux = dst;
  this->do_step(dst, *aux, src, false);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void RelaxationBlockJacobi<MatrixType, InverseNumberType, VectorType>::Tstep
(VectorType       &dst,
 const VectorType &src) const
{
  GrowingVectorMemory<VectorType> mem;
  typename VectorMemory<VectorType>::Pointer aux = mem;
  aux->reinit(dst, false);
  *aux = dst;
  this->do_step(dst, *aux, src, true);
}


//----------------------------------------------------------------------//

template <typename MatrixType, typename InverseNumberType, typename VectorType>
void RelaxationBlockSOR<MatrixType, InverseNumberType, VectorType>::step
(VectorType &dst,
 const VectorType &src) const
{
  this->do_step(dst, dst, src, false);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void RelaxationBlockSOR<MatrixType, InverseNumberType, VectorType>::Tstep
(VectorType       &dst,
 const VectorType &src) const
{
  this->do_step(dst, dst, src, true);
}


//----------------------------------------------------------------------//

template <typename MatrixType, typename InverseNumberType, typename VectorType>
void RelaxationBlockSSOR<MatrixType, InverseNumberType, VectorType>::step
(VectorType       &dst,
 const VectorType &src) const
{
  this->do_step(dst, dst, src, false);
  this->do_step(dst, dst, src, true);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void RelaxationBlockSSOR<MatrixType, InverseNumberType, VectorType>::Tstep
(VectorType       &dst,
 const VectorType &src) const
{
  this->do_step(dst, dst, src, true);
  this->do_step(dst, dst, src, false);
}



DEAL_II_NAMESPACE_CLOSE


#endif
