// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2018 by the deal.II authors
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


#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_block.h>
#include <deal.II/multigrid/mg_transfer_block.templates.h>
#include <deal.II/multigrid/mg_transfer_component.h>
#include <deal.II/multigrid/mg_transfer_component.templates.h>
#include <deal.II/multigrid/multigrid.templates.h>

DEAL_II_NAMESPACE_OPEN


MGTransferBlockBase::MGTransferBlockBase()
  : n_mg_blocks(0)
{}



MGTransferBlockBase::MGTransferBlockBase(const MGConstrainedDoFs &mg_c)
  : n_mg_blocks(0)
  , mg_constrained_dofs(&mg_c)
{}



MGTransferBlockBase::MGTransferBlockBase(
  const AffineConstraints<double> & /*c*/,
  const MGConstrainedDoFs &mg_c)
  : n_mg_blocks(0)
  , mg_constrained_dofs(&mg_c)
{}


template <typename number>
MGTransferBlock<number>::MGTransferBlock()
  : memory(nullptr, typeid(*this).name())
{}


template <typename number>
MGTransferBlock<number>::~MGTransferBlock()
{
  if (memory != nullptr)
    memory = nullptr;
}


template <typename number>
void
MGTransferBlock<number>::initialize(const std::vector<number> &   f,
                                    VectorMemory<Vector<number>> &mem)
{
  factors = f;
  memory  = &mem;
}


template <typename number>
void
MGTransferBlock<number>::prolongate(const unsigned int         to_level,
                                    BlockVector<number> &      dst,
                                    const BlockVector<number> &src) const
{
  Assert((to_level >= 1) && (to_level <= prolongation_matrices.size()),
         ExcIndexRange(to_level, 1, prolongation_matrices.size() + 1));
  Assert(src.n_blocks() == this->n_mg_blocks,
         ExcDimensionMismatch(src.n_blocks(), this->n_mg_blocks));
  Assert(dst.n_blocks() == this->n_mg_blocks,
         ExcDimensionMismatch(dst.n_blocks(), this->n_mg_blocks));

#ifdef DEBUG
  if (this->mg_constrained_dofs != nullptr)
    Assert(this->mg_constrained_dofs->get_user_constraint_matrix(to_level - 1)
               .get_local_lines()
               .size() == 0,
           ExcNotImplemented());
#endif

  // Multiplicate with prolongation
  // matrix, but only those blocks
  // selected.
  for (unsigned int b = 0; b < this->mg_block.size(); ++b)
    {
      if (this->selected[b])
        prolongation_matrices[to_level - 1]->block(b, b).vmult(
          dst.block(this->mg_block[b]), src.block(this->mg_block[b]));
    }
}


template <typename number>
void
MGTransferBlock<number>::restrict_and_add(const unsigned int         from_level,
                                          BlockVector<number> &      dst,
                                          const BlockVector<number> &src) const
{
  Assert((from_level >= 1) && (from_level <= prolongation_matrices.size()),
         ExcIndexRange(from_level, 1, prolongation_matrices.size() + 1));
  Assert(src.n_blocks() == this->n_mg_blocks,
         ExcDimensionMismatch(src.n_blocks(), this->n_mg_blocks));
  Assert(dst.n_blocks() == this->n_mg_blocks,
         ExcDimensionMismatch(dst.n_blocks(), this->n_mg_blocks));

  for (unsigned int b = 0; b < this->mg_block.size(); ++b)
    {
      if (this->selected[b])
        {
          if (factors.size() != 0)
            {
              Assert(memory != nullptr, ExcNotInitialized());
              Vector<number> *aux = memory->alloc();
              aux->reinit(dst.block(this->mg_block[b]));
              prolongation_matrices[from_level - 1]->block(b, b).Tvmult(
                *aux, src.block(this->mg_block[b]));

              dst.block(this->mg_block[b]).add(factors[b], *aux);
              memory->free(aux);
            }
          else
            {
              prolongation_matrices[from_level - 1]->block(b, b).Tvmult_add(
                dst.block(this->mg_block[b]), src.block(this->mg_block[b]));
            }
        }
    }
}



std::size_t
MGTransferComponentBase::memory_consumption() const
{
  std::size_t result = sizeof(*this);
  result += MemoryConsumption::memory_consumption(component_mask) -
            sizeof(ComponentMask);
  result += MemoryConsumption::memory_consumption(target_component) -
            sizeof(mg_target_component);
  result += MemoryConsumption::memory_consumption(sizes) - sizeof(sizes);
  result += MemoryConsumption::memory_consumption(component_start) -
            sizeof(component_start);
  result += MemoryConsumption::memory_consumption(mg_component_start) -
            sizeof(mg_component_start);
  result += MemoryConsumption::memory_consumption(prolongation_sparsities) -
            sizeof(prolongation_sparsities);
  result += MemoryConsumption::memory_consumption(prolongation_matrices) -
            sizeof(prolongation_matrices);
  // TODO:[GK] Add this.
  //   result += MemoryConsumption::memory_consumption(copy_to_and_from_indices)
  //          - sizeof(copy_to_and_from_indices);
  return result;
}


// TODO:[GK] Add all those little vectors.
std::size_t
MGTransferBlockBase::memory_consumption() const
{
  std::size_t result = sizeof(*this);
  result += sizeof(unsigned int) * sizes.size();
  result += MemoryConsumption::memory_consumption(selected) - sizeof(selected);
  result += MemoryConsumption::memory_consumption(mg_block) - sizeof(mg_block);
  result +=
    MemoryConsumption::memory_consumption(block_start) - sizeof(block_start);
  result += MemoryConsumption::memory_consumption(mg_block_start) -
            sizeof(mg_block_start);
  result += MemoryConsumption::memory_consumption(prolongation_sparsities) -
            sizeof(prolongation_sparsities);
  result += MemoryConsumption::memory_consumption(prolongation_matrices) -
            sizeof(prolongation_matrices);
  // TODO:[GK] Add this.
  //   result += MemoryConsumption::memory_consumption(copy_indices)
  //          - sizeof(copy_indices);
  return result;
}


//----------------------------------------------------------------------//

template <typename number>
MGTransferSelect<number>::MGTransferSelect()
  : selected_component(0)
  , mg_selected_component(0)
{}


template <typename number>
MGTransferSelect<number>::MGTransferSelect(const AffineConstraints<double> &c)
  : selected_component(0)
  , mg_selected_component(0)
  , constraints(&c)
{}



template <typename number>
void
MGTransferSelect<number>::prolongate(const unsigned int    to_level,
                                     Vector<number> &      dst,
                                     const Vector<number> &src) const
{
  Assert((to_level >= 1) && (to_level <= prolongation_matrices.size()),
         ExcIndexRange(to_level, 1, prolongation_matrices.size() + 1));

  prolongation_matrices[to_level - 1]
    ->block(mg_target_component[mg_selected_component],
            mg_target_component[mg_selected_component])
    .vmult(dst, src);
}



template <typename number>
void
MGTransferSelect<number>::restrict_and_add(const unsigned int    from_level,
                                           Vector<number> &      dst,
                                           const Vector<number> &src) const
{
  Assert((from_level >= 1) && (from_level <= prolongation_matrices.size()),
         ExcIndexRange(from_level, 1, prolongation_matrices.size() + 1));

  prolongation_matrices[from_level - 1]
    ->block(mg_target_component[mg_selected_component],
            mg_target_component[mg_selected_component])
    .Tvmult_add(dst, src);
}


//----------------------------------------------------------------------//

template <typename number>
MGTransferBlockSelect<number>::MGTransferBlockSelect()
  : selected_block(0)
{}



template <typename number>
MGTransferBlockSelect<number>::MGTransferBlockSelect(
  const MGConstrainedDoFs &mg_c)
  : MGTransferBlockBase(mg_c)
  , selected_block(0)
{}



template <typename number>
MGTransferBlockSelect<number>::MGTransferBlockSelect(
  const AffineConstraints<double> & /*c*/,
  const MGConstrainedDoFs &mg_c)
  : MGTransferBlockBase(mg_c)
  , selected_block(0)
{}



template <typename number>
void
MGTransferBlockSelect<number>::prolongate(const unsigned int    to_level,
                                          Vector<number> &      dst,
                                          const Vector<number> &src) const
{
  Assert((to_level >= 1) && (to_level <= prolongation_matrices.size()),
         ExcIndexRange(to_level, 1, prolongation_matrices.size() + 1));

#ifdef DEBUG
  if (this->mg_constrained_dofs != nullptr)
    Assert(this->mg_constrained_dofs->get_user_constraint_matrix(to_level - 1)
               .get_local_lines()
               .size() == 0,
           ExcNotImplemented());
#endif

  prolongation_matrices[to_level - 1]
    ->block(selected_block, selected_block)
    .vmult(dst, src);
}


template <typename number>
void
MGTransferBlockSelect<number>::restrict_and_add(const unsigned int from_level,
                                                Vector<number> &   dst,
                                                const Vector<number> &src) const
{
  Assert((from_level >= 1) && (from_level <= prolongation_matrices.size()),
         ExcIndexRange(from_level, 1, prolongation_matrices.size() + 1));

  prolongation_matrices[from_level - 1]
    ->block(selected_block, selected_block)
    .Tvmult_add(dst, src);
}



// Explicit instantiations

#include "multigrid.inst"

template class MGTransferBlock<float>;
template class MGTransferBlock<double>;
template class MGTransferSelect<float>;
template class MGTransferSelect<double>;
template class MGTransferBlockSelect<float>;
template class MGTransferBlockSelect<double>;


DEAL_II_NAMESPACE_CLOSE
