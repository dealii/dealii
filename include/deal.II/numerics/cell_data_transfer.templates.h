// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

#ifndef dealii_cell_data_transfer_templates_h
#define dealii_cell_data_transfer_templates_h


#include <deal.II/base/config.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/numerics/cell_data_transfer.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace CellDataTransferImplementation
  {
    template <typename VectorType>
    void
    post_unpack_action(VectorType &out)
    {
      out.compress(::dealii::VectorOperation::insert);
    }

    template <typename value_type>
    void
    post_unpack_action(std::vector<value_type> &)
    {
      // Do nothing for std::vector as VectorType.
    }
  } // namespace CellDataTransferImplementation
} // namespace internal



template <int dim, int spacedim, typename VectorType>
CellDataTransfer<dim, spacedim, VectorType>::CellDataTransfer(
  const Triangulation<dim, spacedim> &               triangulation,
  const std::function<std::vector<value_type>(
    const typename Triangulation<dim, spacedim>::cell_iterator &parent,
    const value_type parent_value)>                  refinement_strategy,
  const std::function<value_type(
    const typename Triangulation<dim, spacedim>::cell_iterator &parent,
    const std::vector<value_type> &children_values)> coarsening_strategy)
  : triangulation(&triangulation, typeid(*this).name())
  , refinement_strategy(refinement_strategy)
  , coarsening_strategy(coarsening_strategy)
  , n_active_cells_pre(numbers::invalid_unsigned_int)
{
  Assert(
    (dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim> *>(
       &triangulation) == nullptr),
    ExcMessage("You are calling the CellDataTransfer class "
               "with a parallel::distributed::Triangulation. "
               "You probably want to use the "
               "parallel::distributed::CellDataTransfer class."));
}



template <int dim, int spacedim, typename VectorType>
void
CellDataTransfer<dim, spacedim, VectorType>::
  prepare_for_coarsening_and_refinement()
{
  // Cleanup previous indices.
  refined_cells_active_index.clear();
  coarsened_cells_active_index.clear();
  persisting_cells_active_index.clear();

  for (const auto &cell : triangulation->active_cell_iterators())
    {
      // Store iterator and active cell index of cells that will be refined.
      if (cell->refine_flag_set())
        {
          refined_cells_active_index.insert({cell, cell->active_cell_index()});
        }
      else if (cell->coarsen_flag_set())
        {
          // Gather the iterator to the parent cell of cells that will be
          // coarsened. Store it together with the active cell indices of all
          // its children.
          Assert(cell->level() > 0, ExcInternalError());
          const auto &parent = cell->parent();

          // Check if the active_cell_indices for the current cell have
          // been determined already.
          if (coarsened_cells_active_index.find(parent) ==
              coarsened_cells_active_index.end())
            {
              std::set<unsigned int> indices_children;
              for (unsigned int child_index = 0;
                   child_index < parent->n_children();
                   ++child_index)
                {
                  const auto &sibling = parent->child(child_index);
                  Assert(sibling->is_active() && sibling->coarsen_flag_set(),
                         typename dealii::Triangulation<
                           dim>::ExcInconsistentCoarseningFlags());

                  indices_children.insert(sibling->active_cell_index());
                }
              AssertDimension(indices_children.size(), parent->n_children());

              coarsened_cells_active_index.insert({parent, indices_children});
            }
        }
      else
        {
          // Store iterator and active cell index of all other cells.
          persisting_cells_active_index.insert(
            {cell, cell->active_cell_index()});
        }
    }

#ifdef DEBUG
  n_active_cells_pre = triangulation->n_active_cells();
#else
  (void)n_active_cells_pre;
#endif
}



template <int dim, int spacedim, typename VectorType>
void
CellDataTransfer<dim, spacedim, VectorType>::unpack(const VectorType &in,
                                                    VectorType &      out)
{
#ifdef DEBUG
  Assert(in.size() == n_active_cells_pre,
         ExcDimensionMismatch(in.size(), n_active_cells_pre));
  Assert(out.size() == triangulation->n_active_cells(),
         ExcDimensionMismatch(out.size(), triangulation->n_active_cells()));
#else
  (void)n_active_cells_pre;
#endif

  // Transfer data of persisting cells.
  for (const auto &persisting : persisting_cells_active_index)
    {
      Assert(persisting.first->is_active(), ExcInternalError());
      out[persisting.first->active_cell_index()] = in[persisting.second];
    }

  // Transfer data of the parent cell to all of its children that it has been
  // refined to.
  std::vector<value_type> children_values;
  for (const auto &refined : refined_cells_active_index)
    {
      // Decide how to handle the previous data.
      children_values = refinement_strategy(refined.first, in[refined.second]);
      Assert(refined.first->n_children() == children_values.size(),
             ExcInternalError());

      // Set values for all children cells.
      for (unsigned int child_index = 0;
           child_index < refined.first->n_children();
           ++child_index)
        {
          const auto &child = refined.first->child(child_index);
          Assert(child->is_active(), ExcInternalError());
          out[child->active_cell_index()] = children_values[child_index];
        }
    }

  // Transfer data from the former children to the cell that they have been
  // coarsened to.
  for (const auto &coarsened : coarsened_cells_active_index)
    {
      // Get previous values of former children.
      children_values.resize(coarsened.second.size());

      auto indices_it = coarsened.second.cbegin();
      auto values_it  = children_values.begin();
      for (; indices_it != coarsened.second.cend(); ++indices_it, ++values_it)
        *values_it = in[*indices_it];
      Assert(values_it == children_values.end(), ExcInternalError());

      // Decide how to handle the previous data.
      const value_type parent_value =
        coarsening_strategy(coarsened.first, children_values);

      // Set value for the parent cell.
      Assert(coarsened.first->is_active(), ExcInternalError());
      out[coarsened.first->active_cell_index()] = parent_value;
    }

  internal::CellDataTransferImplementation::post_unpack_action(out);
}


DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_cell_data_transfer_templates_h */
