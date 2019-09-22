// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_relaxation_block_templates_h
#define dealii_relaxation_block_templates_h

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/relaxation_block.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline RelaxationBlock<MatrixType, InverseNumberType, VectorType>::
  AdditionalData::AdditionalData(
    const double relaxation,
    const bool   invert_diagonal,
    const bool   same_diagonal,
    const typename PreconditionBlockBase<InverseNumberType>::Inversion
                 inversion,
    const double threshold,
    VectorType * temp_ghost_vector)
  : relaxation(relaxation)
  , invert_diagonal(invert_diagonal)
  , same_diagonal(same_diagonal)
  , inversion(inversion)
  , threshold(threshold)
  , temp_ghost_vector(temp_ghost_vector)
{}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline std::size_t
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::AdditionalData::
  memory_consumption() const
{
  std::size_t result =
    sizeof(*this) + block_list.memory_consumption() - sizeof(block_list);
  for (unsigned int i = 0; i < order.size(); ++i)
    result += MemoryConsumption::memory_consumption(order[i]);
  return result;
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline void
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::initialize(
  const MatrixType &    M,
  const AdditionalData &parameters)
{
  Assert(parameters.invert_diagonal, ExcNotImplemented());

  clear();
  //  Assert (M.m() == M.n(), ExcNotQuadratic());
  A               = &M;
  additional_data = &parameters;
  this->inversion = parameters.inversion;

  this->reinit(additional_data->block_list.n_rows(),
               0,
               additional_data->same_diagonal,
               additional_data->inversion);

  if (additional_data->invert_diagonal)
    invert_diagblocks();
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline void
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::clear()
{
  A               = nullptr;
  additional_data = nullptr;
  PreconditionBlockBase<InverseNumberType>::clear();
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline void
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::invert_diagblocks()
{
  if (this->same_diagonal())
    {
      Assert(false, ExcNotImplemented());
    }
  else
    {
      // compute blocks in parallel
      parallel::apply_to_subranges(0,
                                   this->additional_data->block_list.n_rows(),
                                   [this](const size_type block_begin,
                                          const size_type block_end) {
                                     this->block_kernel(block_begin, block_end);
                                   },
                                   16);
    }
  this->inverses_computed(true);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline void
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::block_kernel(
  const size_type block_begin,
  const size_type block_end)
{
  const MatrixType &            M = *(this->A);
  FullMatrix<InverseNumberType> M_cell;

  for (size_type block = block_begin; block < block_end; ++block)
    {
      const size_type bs = this->additional_data->block_list.row_length(block);
      M_cell.reinit(bs, bs);

      // Copy rows for this block into the matrix for the diagonal block
      SparsityPattern::iterator row =
        this->additional_data->block_list.begin(block);
      for (size_type row_cell = 0; row_cell < bs; ++row_cell, ++row)
        {
          for (typename MatrixType::const_iterator entry =
                 M.begin(row->column());
               entry != M.end(row->column());
               ++entry)
            {
              const size_type column = entry->column();
              const size_type col_cell =
                this->additional_data->block_list.row_position(block, column);
              if (col_cell != numbers::invalid_size_type)
                M_cell(row_cell, col_cell) = entry->value();
            }
        }
      // Now M_cell contains the diagonal block. Now store it and its
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
            if (this->additional_data->kernel_size > 0)
              this->inverse_svd(block).compute_inverse_svd_with_kernel(
                this->additional_data->kernel_size);
            else
              this->inverse_svd(block).compute_inverse_svd(
                this->additional_data->threshold);
            break;
          default:
            Assert(false, ExcNotImplemented());
        }
    }
}

namespace internal
{
  /**
   * Default implementation for serial vectors. Here we don't need to make a
   * copy into a ghosted vector, so just return a reference to @p prev.
   */
  template <class VectorType>
  const VectorType &
  prepare_ghost_vector(const VectorType &prev, VectorType *other)
  {
    // If the following Assertion triggers, you either set temp_ghost_vector
    // for a serial computation (don't!), or nobody implemented, instantiated,
    // and tested the parallel version for your vector type.
    Assert(other == nullptr, ExcNotImplemented());
    (void)other;
    return prev;
  }

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * Specialization for Trilinos. Use the ghosted vector.
   */
  template <>
  const TrilinosWrappers::MPI::Vector &
  prepare_ghost_vector(const TrilinosWrappers::MPI::Vector &prev,
                       TrilinosWrappers::MPI::Vector *      other)
  {
    Assert(
      other != nullptr,
      ExcMessage(
        "You need to provide a ghosted vector in RelaxationBlock::AdditionalData::temp_trilinos_ghost_vector."));
    Assert(other->size() == prev.size(), ExcInternalError());

    // import ghost values:
    *other = prev;
    return *other;
  }
#endif // DEAL_II_WITH_TRILINOS
} // end namespace internal

template <typename MatrixType, typename InverseNumberType, typename VectorType>
inline void
RelaxationBlock<MatrixType, InverseNumberType, VectorType>::do_step(
  VectorType &      dst,
  const VectorType &prev,
  const VectorType &src,
  const bool        backward) const
{
  Assert(additional_data->invert_diagonal, ExcNotImplemented());

  const VectorType &ghosted_prev =
    internal::prepare_ghost_vector(prev, additional_data->temp_ghost_vector);

  const MatrixType &                      M = *this->A;
  Vector<typename VectorType::value_type> b_cell, x_cell;

  const bool         permutation_empty = additional_data->order.size() == 0;
  const unsigned int n_permutations =
    (permutation_empty) ? 1U : additional_data->order.size();
  const size_type n_blocks = additional_data->block_list.n_rows();

  if (!permutation_empty)
    for (unsigned int i = 0; i < additional_data->order.size(); ++i)
      AssertDimension(additional_data->order[i].size(), this->size());

  for (unsigned int perm = 0; perm < n_permutations; ++perm)
    {
      for (unsigned int bi = 0; bi < n_blocks; ++bi)
        {
          const unsigned int raw_block = backward ? (n_blocks - bi - 1) : bi;
          const unsigned int block =
            permutation_empty ?
              raw_block :
              (backward ? (additional_data
                             ->order[n_permutations - 1 - perm][raw_block]) :
                          (additional_data->order[perm][raw_block]));

          const size_type bs = additional_data->block_list.row_length(block);

          b_cell.reinit(bs);
          x_cell.reinit(bs);
          // Collect off-diagonal parts
          SparsityPattern::iterator row =
            additional_data->block_list.begin(block);
          for (size_type row_cell = 0; row_cell < bs; ++row_cell, ++row)
            {
              b_cell(row_cell) = src(row->column());
              for (typename MatrixType::const_iterator entry =
                     M.begin(row->column());
                   entry != M.end(row->column());
                   ++entry)
                b_cell(row_cell) -=
                  entry->value() * ghosted_prev(entry->column());
            }
          // Apply inverse diagonal
          this->inverse_vmult(block, x_cell, b_cell);
#ifdef DEBUG
          for (unsigned int i = 0; i < x_cell.size(); ++i)
            {
              AssertIsFinite(x_cell(i));
            }
#endif
          // Store in result vector
          row = additional_data->block_list.begin(block);
          for (size_type row_cell = 0; row_cell < bs; ++row_cell, ++row)
            dst(row->column()) +=
              additional_data->relaxation * x_cell(row_cell);
        }
    }
  dst.compress(dealii::VectorOperation::add);
}


//----------------------------------------------------------------------//

template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockJacobi<MatrixType, InverseNumberType, VectorType>::step(
  VectorType &      dst,
  const VectorType &src) const
{
  GrowingVectorMemory<VectorType>            mem;
  typename VectorMemory<VectorType>::Pointer aux(mem);
  aux->reinit(dst, false);
  *aux = dst;
  this->do_step(dst, *aux, src, false);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockJacobi<MatrixType, InverseNumberType, VectorType>::Tstep(
  VectorType &      dst,
  const VectorType &src) const
{
  GrowingVectorMemory<VectorType>            mem;
  typename VectorMemory<VectorType>::Pointer aux(mem);
  aux->reinit(dst, false);
  *aux = dst;
  this->do_step(dst, *aux, src, true);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockJacobi<MatrixType, InverseNumberType, VectorType>::vmult(
  VectorType &      dst,
  const VectorType &src) const
{
  GrowingVectorMemory<VectorType>            mem;
  typename VectorMemory<VectorType>::Pointer aux(mem);
  dst = 0;
  aux->reinit(dst);
  this->do_step(dst, *aux, src, false);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockJacobi<MatrixType, InverseNumberType, VectorType>::Tvmult(
  VectorType &      dst,
  const VectorType &src) const
{
  GrowingVectorMemory<VectorType>            mem;
  typename VectorMemory<VectorType>::Pointer aux(mem);
  dst = 0;
  aux->reinit(dst);
  this->do_step(dst, *aux, src, true);
}


//----------------------------------------------------------------------//

template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockSOR<MatrixType, InverseNumberType, VectorType>::step(
  VectorType &      dst,
  const VectorType &src) const
{
  this->do_step(dst, dst, src, false);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockSOR<MatrixType, InverseNumberType, VectorType>::Tstep(
  VectorType &      dst,
  const VectorType &src) const
{
  this->do_step(dst, dst, src, true);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockSOR<MatrixType, InverseNumberType, VectorType>::vmult(
  VectorType &      dst,
  const VectorType &src) const
{
  dst = 0;
  this->do_step(dst, dst, src, false);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockSOR<MatrixType, InverseNumberType, VectorType>::Tvmult(
  VectorType &      dst,
  const VectorType &src) const
{
  dst = 0;
  this->do_step(dst, dst, src, true);
}


//----------------------------------------------------------------------//

template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockSSOR<MatrixType, InverseNumberType, VectorType>::step(
  VectorType &      dst,
  const VectorType &src) const
{
  this->do_step(dst, dst, src, false);
  this->do_step(dst, dst, src, true);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockSSOR<MatrixType, InverseNumberType, VectorType>::Tstep(
  VectorType &      dst,
  const VectorType &src) const
{
  this->do_step(dst, dst, src, true);
  this->do_step(dst, dst, src, false);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockSSOR<MatrixType, InverseNumberType, VectorType>::vmult(
  VectorType &      dst,
  const VectorType &src) const
{
  dst = 0;
  this->do_step(dst, dst, src, false);
  this->do_step(dst, dst, src, true);
}


template <typename MatrixType, typename InverseNumberType, typename VectorType>
void
RelaxationBlockSSOR<MatrixType, InverseNumberType, VectorType>::Tvmult(
  VectorType &      dst,
  const VectorType &src) const
{
  dst = 0;
  this->do_step(dst, dst, src, true);
  this->do_step(dst, dst, src, false);
}



DEAL_II_NAMESPACE_CLOSE


#endif
