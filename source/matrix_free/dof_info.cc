// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/matrix_free/dof_info.templates.h>
#include <deal.II/matrix_free/vector_data_exchange.h>

#include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    // ensure that the type defined in both dof_info.h and
    // hanging_nodes_internal.h is consistent
    static_assert(std::is_same_v<compressed_constraint_kind, std::uint8_t>,
                  "Unexpected type for compressed hanging node indicators!");



    DoFInfo::DoFInfo()
    {
      clear();
    }



    void
    DoFInfo::clear()
    {
      row_starts.clear();
      dof_indices.clear();
      constraint_indicator.clear();
      vector_partitioner.reset();
      ghost_dofs.clear();
      dofs_per_cell.clear();
      dofs_per_face.clear();
      vectorization_length       = 1;
      global_base_element_offset = 0;
      n_base_elements            = 0;
      n_components.clear();
      start_components.clear();
      row_starts_plain_indices.clear();
      plain_dof_indices.clear();
      dof_indices_interleaved.clear();
      for (unsigned int i = 0; i < 3; ++i)
        {
          index_storage_variants[i].clear();
          dof_indices_contiguous[i].clear();
          dof_indices_interleave_strides[i].clear();
          n_vectorization_lanes_filled[i].clear();
        }
      store_plain_indices = false;
      cell_active_fe_index.clear();
      max_fe_index = 0;
      fe_index_conversion.clear();
    }



    void
    DoFInfo::get_dof_indices_on_cell_batch(std::vector<unsigned int> &my_rows,
                                           const unsigned int         cell,
                                           const bool apply_constraints) const
    {
      const unsigned int n_fe_components = start_components.back();
      const unsigned int fe_index =
        dofs_per_cell.size() == 1 ? 0 : cell_active_fe_index[cell];
      const unsigned int dofs_this_cell = dofs_per_cell[fe_index];

      const unsigned int n_vectorization  = vectorization_length;
      constexpr auto     dof_access_index = dof_access_cell;
      AssertIndexRange(cell,
                       n_vectorization_lanes_filled[dof_access_index].size());
      const unsigned int n_vectorization_actual =
        n_vectorization_lanes_filled[dof_access_index][cell];

      // we might have constraints, so the final number
      // of indices is not known a priori.
      // conservatively reserve the maximum without constraints
      my_rows.reserve(n_vectorization * dofs_this_cell);
      my_rows.resize(0);
      unsigned int total_size = 0;
      for (unsigned int v = 0; v < n_vectorization_actual; ++v)
        {
          const unsigned int ib =
            (cell * n_vectorization + v) * n_fe_components;
          const unsigned int ie =
            (cell * n_vectorization + v + 1) * n_fe_components;

          // figure out constraints by comparing constraint_indicator row
          // shift for this cell within the block as compared to the next
          // one
          const bool has_constraints =
            (hanging_node_constraint_masks.size() != 0 &&
             hanging_node_constraint_masks_comp.size() != 0 &&
             hanging_node_constraint_masks[cell * n_vectorization + v] !=
               unconstrained_compressed_constraint_kind &&
             hanging_node_constraint_masks_comp[fe_index][0 /*TODO*/]) ||
            (row_starts[ib].second != row_starts[ib + n_fe_components].second);

          auto do_copy = [&](const unsigned int *begin,
                             const unsigned int *end) {
            const unsigned int shift = total_size;
            total_size += (end - begin);
            my_rows.resize(total_size);
            std::copy(begin, end, my_rows.begin() + shift);
          };

          if (!has_constraints || apply_constraints)
            {
              const unsigned int *begin =
                dof_indices.data() + row_starts[ib].first;
              const unsigned int *end =
                dof_indices.data() + row_starts[ie].first;
              do_copy(begin, end);
            }
          else
            {
              Assert(row_starts_plain_indices[cell * n_vectorization + v] !=
                       numbers::invalid_unsigned_int,
                     ExcNotInitialized());
              const unsigned int *begin =
                plain_dof_indices.data() +
                row_starts_plain_indices[cell * n_vectorization + v];
              const unsigned int *end = begin + dofs_this_cell;
              do_copy(begin, end);
            }
        }
    }



    void
    DoFInfo::assign_ghosts(const std::vector<unsigned int> &boundary_cells,
                           const MPI_Comm                   communicator_sm,
                           const bool use_vector_data_exchanger_full)
    {
      Assert(boundary_cells.size() < row_starts.size(), ExcInternalError());

      // sort ghost dofs and compress out duplicates
      const unsigned int n_owned  = (vector_partitioner->local_range().second -
                                    vector_partitioner->local_range().first);
      const std::size_t  n_ghosts = ghost_dofs.size();
      if constexpr (running_in_debug_mode())
        {
          for (const auto dof_index : dof_indices)
            AssertIndexRange(dof_index, n_owned + n_ghosts);
        }

      const unsigned int        n_components = start_components.back();
      std::vector<unsigned int> ghost_numbering(n_ghosts);
      IndexSet                  ghost_indices(vector_partitioner->size());
      if (n_ghosts > 0)
        {
          unsigned int n_unique_ghosts = 0;
          // since we need to go back to the local_to_global indices and
          // replace the temporary numbering of ghosts by the real number in
          // the index set, we need to store these values
          std::vector<std::pair<types::global_dof_index, unsigned int>>
            ghost_origin(n_ghosts);
          for (std::size_t i = 0; i < n_ghosts; ++i)
            {
              ghost_origin[i].first  = ghost_dofs[i];
              ghost_origin[i].second = i;
            }
          std::sort(ghost_origin.begin(), ghost_origin.end());

          types::global_dof_index last_contiguous_start = ghost_origin[0].first;
          ghost_numbering[ghost_origin[0].second]       = 0;
          for (std::size_t i = 1; i < n_ghosts; ++i)
            {
              if (ghost_origin[i].first > ghost_origin[i - 1].first + 1)
                {
                  ghost_indices.add_range(last_contiguous_start,
                                          ghost_origin[i - 1].first + 1);
                  last_contiguous_start = ghost_origin[i].first;
                }
              if (ghost_origin[i].first > ghost_origin[i - 1].first)
                ++n_unique_ghosts;
              ghost_numbering[ghost_origin[i].second] = n_unique_ghosts;
            }
          ++n_unique_ghosts;
          ghost_indices.add_range(last_contiguous_start,
                                  ghost_origin.back().first + 1);
          ghost_indices.compress();

          // make sure that we got the correct local numbering of the ghost
          // dofs. the ghost index set should store the same number
          {
            AssertDimension(n_unique_ghosts, ghost_indices.n_elements());
            for (std::size_t i = 0; i < n_ghosts; ++i)
              Assert(ghost_numbering[i] ==
                       ghost_indices.index_within_set(ghost_dofs[i]),
                     ExcInternalError());
          }

          // apply correct numbering for ghost indices: We previously just
          // enumerated them according to their appearance in the
          // local_to_global structure. Above, we derived a relation between
          // this enumeration and the actual number
          const unsigned int n_boundary_cells = boundary_cells.size();
          for (unsigned int i = 0; i < n_boundary_cells; ++i)
            {
              unsigned int *data_ptr =
                dof_indices.data() +
                row_starts[boundary_cells[i] * n_components].first;
              const unsigned int *row_end =
                dof_indices.data() +
                row_starts[(boundary_cells[i] + 1) * n_components].first;
              for (; data_ptr != row_end; ++data_ptr)
                *data_ptr = ((*data_ptr < n_owned) ?
                               *data_ptr :
                               n_owned + ghost_numbering[*data_ptr - n_owned]);

              // now the same procedure for plain indices
              if (store_plain_indices == true)
                {
                  bool has_hanging_nodes = false;

                  const unsigned int fe_index =
                    (cell_active_fe_index.empty() ||
                     dofs_per_cell.size() == 1) ?
                      0 :
                      cell_active_fe_index[boundary_cells[i]];
                  AssertIndexRange(fe_index, dofs_per_cell.size());

                  if (hanging_node_constraint_masks.size() > 0 &&
                      hanging_node_constraint_masks_comp.size() > 0 &&
                      hanging_node_constraint_masks[boundary_cells[i]] !=
                        unconstrained_compressed_constraint_kind)
                    for (unsigned int comp = 0; comp < n_components; ++comp)
                      has_hanging_nodes |=
                        hanging_node_constraint_masks_comp[fe_index][comp];

                  if (has_hanging_nodes ||
                      row_starts[boundary_cells[i] * n_components].second !=
                        row_starts[(boundary_cells[i] + 1) * n_components]
                          .second)
                    {
                      unsigned int *data_ptr =
                        plain_dof_indices.data() +
                        row_starts_plain_indices[boundary_cells[i]];
                      const unsigned int *row_end =
                        data_ptr + dofs_per_cell[fe_index];
                      for (; data_ptr != row_end; ++data_ptr)
                        *data_ptr =
                          ((*data_ptr < n_owned) ?
                             *data_ptr :
                             n_owned + ghost_numbering[*data_ptr - n_owned]);
                    }
                }
            }
        }

      std::vector<types::global_dof_index> empty;
      ghost_dofs.swap(empty);

      // set the ghost indices now. need to cast away constness here, but that
      // is uncritical since we reset the Partitioner in the same initialize
      // call as this call here.
      Utilities::MPI::Partitioner *vec_part =
        const_cast<Utilities::MPI::Partitioner *>(vector_partitioner.get());
      vec_part->set_ghost_indices(ghost_indices);

      if (use_vector_data_exchanger_full == false)
        vector_exchanger = std::make_shared<
          MatrixFreeFunctions::VectorDataExchange::PartitionerWrapper>(
          vector_partitioner);
      else
        vector_exchanger =
          std::make_shared<MatrixFreeFunctions::VectorDataExchange::Full>(
            vector_partitioner, communicator_sm);
    }



    void
    DoFInfo::reorder_cells(
      const TaskInfo                   &task_info,
      const std::vector<unsigned int>  &renumbering,
      const std::vector<unsigned int>  &constraint_pool_row_index,
      const std::vector<unsigned char> &irregular_cells)
    {
      // first reorder the active FE index.
      const bool have_hp = dofs_per_cell.size() > 1;
      if (cell_active_fe_index.size() > 0)
        {
          std::vector<unsigned int> new_active_fe_index;
          new_active_fe_index.reserve(task_info.cell_partition_data.back());
          unsigned int position_cell = 0;
          for (unsigned int cell = 0;
               cell < task_info.cell_partition_data.back();
               ++cell)
            {
              const unsigned int n_comp =
                (irregular_cells[cell] > 0 ? irregular_cells[cell] :
                                             vectorization_length);

              // take maximum FE index among the ones present (we might have
              // lumped some lower indices into higher ones)
              unsigned int fe_index =
                cell_active_fe_index[renumbering[position_cell]];
              for (unsigned int j = 1; j < n_comp; ++j)
                fe_index = std::max(
                  fe_index,
                  cell_active_fe_index[renumbering[position_cell + j]]);

              new_active_fe_index.push_back(fe_index);
              position_cell += n_comp;
            }
          std::swap(new_active_fe_index, cell_active_fe_index);
        }
      if (have_hp)
        AssertDimension(cell_active_fe_index.size(),
                        task_info.cell_partition_data.back());

      const unsigned int n_components = start_components.back();

      std::vector<std::pair<unsigned int, unsigned int>> new_row_starts(
        vectorization_length * n_components *
          task_info.cell_partition_data.back() +
        1);
      std::vector<unsigned int> new_dof_indices;
      std::vector<std::pair<unsigned short, unsigned short>>
                                new_constraint_indicator;
      std::vector<unsigned int> new_plain_indices, new_rowstart_plain;
      unsigned int              position_cell = 0;
      new_dof_indices.reserve(dof_indices.size());
      new_constraint_indicator.reserve(constraint_indicator.size());

      std::vector<compressed_constraint_kind> new_hanging_node_constraint_masks;
      new_hanging_node_constraint_masks.reserve(
        hanging_node_constraint_masks.size());

      if (store_plain_indices == true)
        {
          new_rowstart_plain.resize(vectorization_length *
                                        task_info.cell_partition_data.back() +
                                      1,
                                    numbers::invalid_unsigned_int);
          new_plain_indices.reserve(plain_dof_indices.size());
        }

      // copy the indices and the constraint indicators to the new data field,
      // where we will go through the cells in the renumbered way. in case the
      // vectorization length does not exactly match up, we fill invalid
      // numbers to the rowstart data. for contiguous cell indices, we skip
      // the rowstarts field completely and directly go into the
      // new_dof_indices field (this layout is used in FEEvaluation).
      for (unsigned int i = 0; i < task_info.cell_partition_data.back(); ++i)
        {
          const unsigned int n_vect =
            (irregular_cells[i] > 0 ? irregular_cells[i] :
                                      vectorization_length);
          const unsigned int dofs_per_cell =
            have_hp ? this->dofs_per_cell[cell_active_fe_index[i]] :
                      this->dofs_per_cell[0];

          for (unsigned int j = 0; j < n_vect; ++j)
            {
              const unsigned int cell_no =
                renumbering[position_cell + j] * n_components;

              bool has_hanging_nodes = false;

              if (hanging_node_constraint_masks.size() > 0 &&
                  hanging_node_constraint_masks_comp.size() > 0)
                {
                  const auto mask =
                    hanging_node_constraint_masks[renumbering[position_cell +
                                                              j]];
                  new_hanging_node_constraint_masks.push_back(mask);

                  if (mask != unconstrained_compressed_constraint_kind)
                    for (unsigned int comp = 0; comp < n_components; ++comp)
                      has_hanging_nodes |= hanging_node_constraint_masks_comp
                        [have_hp ? cell_active_fe_index[i] : 0][comp];
                }

              for (unsigned int comp = 0; comp < n_components; ++comp)
                {
                  new_row_starts[(i * vectorization_length + j) * n_components +
                                 comp]
                    .first = new_dof_indices.size();
                  new_row_starts[(i * vectorization_length + j) * n_components +
                                 comp]
                    .second = new_constraint_indicator.size();

                  new_dof_indices.insert(
                    new_dof_indices.end(),
                    dof_indices.data() + row_starts[cell_no + comp].first,
                    dof_indices.data() + row_starts[cell_no + comp + 1].first);
                  for (unsigned int index = row_starts[cell_no + comp].second;
                       index != row_starts[cell_no + comp + 1].second;
                       ++index)
                    new_constraint_indicator.push_back(
                      constraint_indicator[index]);
                }
              if (store_plain_indices &&
                  ((row_starts[cell_no].second !=
                    row_starts[cell_no + n_components].second) ||
                   has_hanging_nodes))
                {
                  new_rowstart_plain[i * vectorization_length + j] =
                    new_plain_indices.size();
                  new_plain_indices.insert(
                    new_plain_indices.end(),
                    plain_dof_indices.data() +
                      row_starts_plain_indices[cell_no / n_components],
                    plain_dof_indices.data() +
                      row_starts_plain_indices[cell_no / n_components] +
                      dofs_per_cell);
                }
            }
          for (unsigned int j = n_vect; j < vectorization_length; ++j)
            for (unsigned int comp = 0; comp < n_components; ++comp)
              {
                new_row_starts[(i * vectorization_length + j) * n_components +
                               comp]
                  .first = new_dof_indices.size();
                new_row_starts[(i * vectorization_length + j) * n_components +
                               comp]
                  .second = new_constraint_indicator.size();
              }

          for (unsigned int j = n_vect; j < vectorization_length; ++j)
            if (hanging_node_constraint_masks.size() > 0)
              new_hanging_node_constraint_masks.push_back(
                unconstrained_compressed_constraint_kind);

          position_cell += n_vect;
        }
      AssertDimension(position_cell * n_components + 1, row_starts.size());

      AssertDimension(dof_indices.size(), new_dof_indices.size());
      new_row_starts[task_info.cell_partition_data.back() *
                     vectorization_length * n_components]
        .first = new_dof_indices.size();
      new_row_starts[task_info.cell_partition_data.back() *
                     vectorization_length * n_components]
        .second = new_constraint_indicator.size();

      AssertDimension(constraint_indicator.size(),
                      new_constraint_indicator.size());

      new_row_starts.swap(row_starts);
      new_dof_indices.swap(dof_indices);
      new_constraint_indicator.swap(constraint_indicator);
      new_plain_indices.swap(plain_dof_indices);
      new_rowstart_plain.swap(row_starts_plain_indices);
      new_hanging_node_constraint_masks.swap(hanging_node_constraint_masks);

      if constexpr (running_in_debug_mode())
        {
          // sanity check 1: all indices should be smaller than the number of
          // dofs locally owned plus the number of ghosts
          const unsigned int index_range =
            (vector_partitioner->local_range().second -
             vector_partitioner->local_range().first) +
            vector_partitioner->ghost_indices().n_elements();
          for (const auto dof_index : dof_indices)
            AssertIndexRange(dof_index, index_range);

          // sanity check 2: for the constraint indicators, the first index
          // should be smaller than the number of indices in the row, and the
          // second index should be smaller than the number of constraints in
          // the constraint pool.
          for (unsigned int row = 0; row < task_info.cell_partition_data.back();
               ++row)
            {
              const unsigned int row_length_ind =
                row_starts[(row * vectorization_length + 1) * n_components]
                  .first -
                row_starts[row * vectorization_length * n_components].first;
              AssertIndexRange(
                row_starts[(row * vectorization_length + 1) * n_components]
                  .second,
                constraint_indicator.size() + 1);
              const std::pair<unsigned short, unsigned short>
                *con_it =
                  constraint_indicator.data() +
                  row_starts[row * vectorization_length * n_components].second,
                *end_con =
                  constraint_indicator.data() +
                  row_starts[(row * vectorization_length + 1) * n_components]
                    .second;
              for (; con_it != end_con; ++con_it)
                {
                  AssertIndexRange(con_it->first, row_length_ind + 1);
                  AssertIndexRange(con_it->second,
                                   constraint_pool_row_index.size() - 1);
                }
            }

          // sanity check 3: check the number of cells once again
          unsigned int n_active_cells = 0;
          for (unsigned int c = 0;
               c < *(task_info.cell_partition_data.end() - 2);
               ++c)
            if (irregular_cells[c] > 0)
              n_active_cells += irregular_cells[c];
            else
              n_active_cells += vectorization_length;
          AssertDimension(n_active_cells, task_info.n_active_cells);
        }

      compute_cell_index_compression(irregular_cells);
    }



    void
    DoFInfo::compute_cell_index_compression(
      const std::vector<unsigned char> &irregular_cells)
    {
      const bool         have_hp      = dofs_per_cell.size() > 1;
      const unsigned int n_components = start_components.back();

      Assert(vectorization_length == 1 ||
               row_starts.size() % vectorization_length == 1,
             ExcInternalError());
      if (vectorization_length > 1)
        AssertDimension(row_starts.size() / vectorization_length / n_components,
                        irregular_cells.size());
      index_storage_variants[dof_access_cell].resize(
        irregular_cells.size(), IndexStorageVariants::full);
      n_vectorization_lanes_filled[dof_access_cell].resize(
        irregular_cells.size());
      for (unsigned int i = 0; i < irregular_cells.size(); ++i)
        if (irregular_cells[i] > 0)
          n_vectorization_lanes_filled[dof_access_cell][i] = irregular_cells[i];
        else
          n_vectorization_lanes_filled[dof_access_cell][i] =
            vectorization_length;

      dof_indices_contiguous[dof_access_cell].resize(
        irregular_cells.size() * vectorization_length,
        numbers::invalid_unsigned_int);
      dof_indices_interleaved.resize(dof_indices.size(),
                                     numbers::invalid_unsigned_int);
      dof_indices_interleave_strides[dof_access_cell].resize(
        irregular_cells.size() * vectorization_length,
        numbers::invalid_unsigned_int);

      std::vector<unsigned int> index_kinds(
        static_cast<unsigned int>(
          IndexStorageVariants::interleaved_contiguous_mixed_strides) +
        1);
      std::vector<unsigned int> offsets(vectorization_length);
      for (unsigned int i = 0; i < irregular_cells.size(); ++i)
        {
          const unsigned int ndofs =
            dofs_per_cell[have_hp ? cell_active_fe_index[i] : 0];
          const unsigned int n_comp =
            n_vectorization_lanes_filled[dof_access_cell][i];

          // check 1: Check if there are constraints -> no compression possible
          bool has_constraints = false;
          for (unsigned int j = 0; j < n_comp; ++j)
            {
              const unsigned int cell_no = i * vectorization_length + j;
              if (row_starts[cell_no * n_components].second !=
                  row_starts[(cell_no + 1) * n_components].second)
                {
                  has_constraints = true;
                  break;
                }
            }
          if (has_constraints)
            index_storage_variants[dof_access_cell][i] =
              IndexStorageVariants::full;
          else
            {
              bool indices_are_contiguous = (ndofs > 0);
              for (unsigned int j = 0; j < n_comp; ++j)
                {
                  const unsigned int  cell_no = i * vectorization_length + j;
                  const unsigned int *dof_indices =
                    this->dof_indices.data() +
                    row_starts[cell_no * n_components].first;
                  AssertDimension(
                    ndofs,
                    row_starts[(cell_no + 1) * n_components].first -
                      row_starts[cell_no * n_components].first);
                  for (unsigned int i = 1; i < ndofs; ++i)
                    if (dof_indices[i] != dof_indices[0] + i)
                      {
                        indices_are_contiguous = false;
                        break;
                      }
                }

              bool indices_are_interleaved_and_contiguous =
                (ndofs > 1 && n_comp == vectorization_length);

              {
                const unsigned int *dof_indices =
                  this->dof_indices.data() +
                  row_starts[i * vectorization_length * n_components].first;
                for (unsigned int k = 0;
                     k < ndofs && indices_are_interleaved_and_contiguous;
                     ++k)
                  for (unsigned int j = 0; j < n_comp; ++j)
                    if (dof_indices[j * ndofs + k] !=
                        dof_indices[0] + k * n_comp + j)
                      {
                        indices_are_interleaved_and_contiguous = false;
                        break;
                      }
              }

              if (indices_are_contiguous ||
                  indices_are_interleaved_and_contiguous)
                {
                  for (unsigned int j = 0; j < n_comp; ++j)
                    {
                      const unsigned int start_index =
                        row_starts[(i * vectorization_length + j) *
                                   n_components]
                          .first;
                      AssertIndexRange(start_index, dof_indices.size());
                      dof_indices_contiguous[dof_access_cell]
                                            [i * vectorization_length + j] =
                                              this->dof_indices.empty() ?
                                                0 :
                                                this->dof_indices[start_index];
                    }
                }

              if (indices_are_interleaved_and_contiguous)
                {
                  Assert(n_comp == vectorization_length, ExcInternalError());
                  index_storage_variants[dof_access_cell][i] =
                    IndexStorageVariants::interleaved_contiguous;
                  for (unsigned int j = 0; j < n_comp; ++j)
                    dof_indices_interleave_strides[2][i * vectorization_length +
                                                      j] = n_comp;
                }
              else if (indices_are_contiguous)
                {
                  index_storage_variants[dof_access_cell][i] =
                    IndexStorageVariants::contiguous;
                  for (unsigned int j = 0; j < n_comp; ++j)
                    dof_indices_interleave_strides[2][i * vectorization_length +
                                                      j] = 1;
                }
              else if (ndofs > 0)
                {
                  int                 indices_are_interleaved_and_mixed = 2;
                  const unsigned int *dof_indices =
                    &this->dof_indices[row_starts[i * vectorization_length *
                                                  n_components]
                                         .first];
                  for (unsigned int j = 0; j < n_comp; ++j)
                    offsets[j] =
                      dof_indices[j * ndofs + 1] - dof_indices[j * ndofs];
                  for (unsigned int k = 0;
                       k < ndofs && indices_are_interleaved_and_mixed != 0;
                       ++k)
                    for (unsigned int j = 0; j < n_comp; ++j)
                      // the first if case is to avoid negative offsets
                      // (invalid)
                      if (dof_indices[j * ndofs + 1] < dof_indices[j * ndofs] ||
                          dof_indices[j * ndofs + k] !=
                            dof_indices[j * ndofs] + k * offsets[j])
                        {
                          indices_are_interleaved_and_mixed = 0;
                          break;
                        }
                  if (indices_are_interleaved_and_mixed == 2)
                    {
                      for (unsigned int j = 0; j < n_comp; ++j)
                        dof_indices_interleave_strides
                          [dof_access_cell][i * vectorization_length + j] =
                            offsets[j];
                      for (unsigned int j = 0; j < n_comp; ++j)
                        dof_indices_contiguous[dof_access_cell]
                                              [i * vectorization_length + j] =
                                                dof_indices[j * ndofs];
                      for (unsigned int j = 0; j < n_comp; ++j)
                        if (offsets[j] != vectorization_length)
                          {
                            indices_are_interleaved_and_mixed = 1;
                            break;
                          }
                      if (indices_are_interleaved_and_mixed == 1 ||
                          n_comp != vectorization_length)
                        index_storage_variants[dof_access_cell][i] =
                          IndexStorageVariants::
                            interleaved_contiguous_mixed_strides;
                      else
                        index_storage_variants[dof_access_cell][i] =
                          IndexStorageVariants::interleaved_contiguous_strided;
                    }
                  else
                    {
                      const unsigned int *dof_indices =
                        this->dof_indices.data() +
                        row_starts[i * vectorization_length * n_components]
                          .first;
                      if (n_comp == vectorization_length)
                        index_storage_variants[dof_access_cell][i] =
                          IndexStorageVariants::interleaved;
                      else
                        index_storage_variants[dof_access_cell][i] =
                          IndexStorageVariants::full;

                      // do not use interleaved storage if two entries within
                      // vectorized array point to the same index (scatter not
                      // possible due to race condition)
                      for (const unsigned int *indices = dof_indices;
                           indices != dof_indices + ndofs;
                           ++indices)
                        {
                          bool         is_sorted = true;
                          unsigned int previous  = indices[0];
                          for (unsigned int l = 1; l < n_comp; ++l)
                            {
                              const unsigned int current = indices[l * ndofs];
                              if (current <= previous)
                                is_sorted = false;

                              // the simple check failed, must compare all
                              // indices manually - due to short sizes this
                              // O(n^2) algorithm is better than sorting
                              if (!is_sorted)
                                for (unsigned int j = 0; j < l; ++j)
                                  if (indices[j * ndofs] == current)
                                    {
                                      index_storage_variants
                                        [dof_access_cell][i] =
                                          IndexStorageVariants::full;
                                      break;
                                    }
                              previous = current;
                            }
                        }
                    }
                }
              else // ndofs == 0
                index_storage_variants[dof_access_cell][i] =
                  IndexStorageVariants::full;
            }
          index_kinds[static_cast<unsigned int>(
            index_storage_variants[dof_access_cell][i])]++;
        }

      // Cleanup phase: we want to avoid single cells with different properties
      // than the bulk of the domain in order to avoid extra checks in the face
      // identification.

      // Step 1: check whether the interleaved indices were only assigned to
      // the single cell within a vectorized array.
      auto fix_single_interleaved_indices =
        [&](const IndexStorageVariants variant) {
          if (index_kinds[static_cast<unsigned int>(
                IndexStorageVariants::interleaved_contiguous_mixed_strides)] >
                0 &&
              index_kinds[static_cast<unsigned int>(variant)] > 0)
            for (unsigned int i = 0; i < irregular_cells.size(); ++i)
              {
                if (index_storage_variants[dof_access_cell][i] ==
                      IndexStorageVariants::
                        interleaved_contiguous_mixed_strides &&
                    n_vectorization_lanes_filled[dof_access_cell][i] == 1 &&
                    (variant != IndexStorageVariants::contiguous ||
                     dof_indices_interleave_strides[dof_access_cell]
                                                   [i * vectorization_length] ==
                       1))
                  {
                    index_storage_variants[dof_access_cell][i] = variant;
                    index_kinds[static_cast<unsigned int>(
                      IndexStorageVariants::
                        interleaved_contiguous_mixed_strides)]--;
                    index_kinds[static_cast<unsigned int>(variant)]++;
                  }
              }
        };

      fix_single_interleaved_indices(IndexStorageVariants::full);
      fix_single_interleaved_indices(IndexStorageVariants::contiguous);
      fix_single_interleaved_indices(IndexStorageVariants::interleaved);

      unsigned int n_interleaved =
        index_kinds[static_cast<unsigned int>(
          IndexStorageVariants::interleaved_contiguous)] +
        index_kinds[static_cast<unsigned int>(
          IndexStorageVariants::interleaved_contiguous_strided)] +
        index_kinds[static_cast<unsigned int>(
          IndexStorageVariants::interleaved_contiguous_mixed_strides)];

      // Step 2: fix single contiguous cell among others with interleaved
      // storage
      if (n_interleaved > 0 && index_kinds[static_cast<unsigned int>(
                                 IndexStorageVariants::contiguous)] > 0)
        for (unsigned int i = 0; i < irregular_cells.size(); ++i)
          if (index_storage_variants[dof_access_cell][i] ==
              IndexStorageVariants::contiguous)
            {
              index_storage_variants[dof_access_cell][i] =
                IndexStorageVariants::interleaved_contiguous_mixed_strides;
              index_kinds[static_cast<unsigned int>(
                IndexStorageVariants::contiguous)]--;
              index_kinds[static_cast<unsigned int>(
                IndexStorageVariants::interleaved_contiguous_mixed_strides)]++;
            }

      // Step 3: Interleaved cells are left but also some non-contiguous ones
      // -> revert all to full storage
      if (n_interleaved > 0 &&
          index_kinds[static_cast<unsigned int>(IndexStorageVariants::full)] +
              index_kinds[static_cast<unsigned int>(
                IndexStorageVariants::interleaved)] >
            0)
        for (unsigned int i = 0; i < irregular_cells.size(); ++i)
          if (index_storage_variants[dof_access_cell][i] >
              IndexStorageVariants::contiguous)
            {
              index_kinds[static_cast<unsigned int>(
                index_storage_variants[2][i])]--;
              if (n_vectorization_lanes_filled[dof_access_cell][i] ==
                  vectorization_length)
                index_storage_variants[dof_access_cell][i] =
                  IndexStorageVariants::interleaved;
              else
                index_storage_variants[dof_access_cell][i] =
                  IndexStorageVariants::full;
              index_kinds[static_cast<unsigned int>(
                index_storage_variants[dof_access_cell][i])]++;
            }

      // Step 4: Copy the interleaved indices into their own data structure
      for (unsigned int i = 0; i < irregular_cells.size(); ++i)
        if (index_storage_variants[dof_access_cell][i] ==
            IndexStorageVariants::interleaved)
          {
            if (n_vectorization_lanes_filled[dof_access_cell][i] <
                vectorization_length)
              {
                index_storage_variants[dof_access_cell][i] =
                  IndexStorageVariants::full;
                continue;
              }
            const unsigned int ndofs =
              dofs_per_cell[have_hp ? cell_active_fe_index[i] : 0];
            const unsigned int *dof_indices =
              &this->dof_indices
                 [row_starts[i * vectorization_length * n_components].first];
            unsigned int *interleaved_dof_indices =
              &this->dof_indices_interleaved
                 [row_starts[i * vectorization_length * n_components].first];
            AssertDimension(this->dof_indices.size(),
                            this->dof_indices_interleaved.size());
            AssertDimension(n_vectorization_lanes_filled[dof_access_cell][i],
                            vectorization_length);
            AssertIndexRange(
              row_starts[i * vectorization_length * n_components].first,
              this->dof_indices_interleaved.size() + 1);
            AssertIndexRange(
              row_starts[i * vectorization_length * n_components].first +
                ndofs * vectorization_length,
              this->dof_indices_interleaved.size() + 1);
            for (unsigned int k = 0; k < ndofs; ++k)
              {
                const unsigned int *my_dof_indices = dof_indices + k;
                const unsigned int *end =
                  interleaved_dof_indices + vectorization_length;
                for (; interleaved_dof_indices != end;
                     ++interleaved_dof_indices, my_dof_indices += ndofs)
                  *interleaved_dof_indices = *my_dof_indices;
              }
          }
    }



    void
    DoFInfo::compute_tight_partitioners(
      const Table<2, ShapeInfo<double>>        &shape_info,
      const unsigned int                        n_owned_cells,
      const unsigned int                        n_lanes,
      const std::vector<FaceToCellTopology<1>> &inner_faces,
      const std::vector<FaceToCellTopology<1>> &ghosted_faces,
      const bool                                fill_cell_centric,
      const MPI_Comm                            communicator_sm,
      const bool                                use_vector_data_exchanger_full)
    {
      const Utilities::MPI::Partitioner &part = *vector_partitioner;

      // partitioner 0: no face integrals, simply use the indices present
      // on the cells
      std::vector<types::global_dof_index> ghost_indices;
      {
        const unsigned int n_components = start_components.back();
        for (unsigned int cell = 0; cell < n_owned_cells; ++cell)
          {
            for (unsigned int i = row_starts[cell * n_components].first;
                 i < row_starts[(cell + 1) * n_components].first;
                 ++i)
              if (dof_indices[i] >= part.locally_owned_size())
                ghost_indices.push_back(part.local_to_global(dof_indices[i]));

            const unsigned int fe_index =
              dofs_per_cell.size() == 1 ? 0 :
                                          cell_active_fe_index[cell / n_lanes];
            const unsigned int dofs_this_cell = dofs_per_cell[fe_index];

            for (unsigned int i = row_starts_plain_indices[cell];
                 i < row_starts_plain_indices[cell] + dofs_this_cell;
                 ++i)
              if (plain_dof_indices[i] >= part.locally_owned_size())
                ghost_indices.push_back(
                  part.local_to_global(plain_dof_indices[i]));
          }
        std::sort(ghost_indices.begin(), ghost_indices.end());
        IndexSet compressed_set(part.size());
        compressed_set.add_indices(ghost_indices.begin(), ghost_indices.end());
        compressed_set.subtract_set(part.locally_owned_range());
        const bool all_ghosts_equal =
          Utilities::MPI::min<int>(static_cast<int>(
                                     compressed_set.n_elements() ==
                                     part.ghost_indices().n_elements()),
                                   part.get_mpi_communicator()) != 0;

        std::shared_ptr<const Utilities::MPI::Partitioner> temp_0;

        if (all_ghosts_equal)
          temp_0 = vector_partitioner;
        else
          {
            temp_0 = std::make_shared<Utilities::MPI::Partitioner>(
              part.locally_owned_range(), part.get_mpi_communicator());
            const_cast<Utilities::MPI::Partitioner *>(temp_0.get())
              ->set_ghost_indices(compressed_set, part.ghost_indices());
          }

        if (use_vector_data_exchanger_full == false)
          vector_exchanger_face_variants[0] = std::make_shared<
            MatrixFreeFunctions::VectorDataExchange::PartitionerWrapper>(
            temp_0);
        else
          vector_exchanger_face_variants[0] =
            std::make_shared<MatrixFreeFunctions::VectorDataExchange::Full>(
              temp_0, communicator_sm);
      }

      // construct a numbering of faces
      std::vector<FaceToCellTopology<1>> all_faces(inner_faces);
      all_faces.insert(all_faces.end(),
                       ghosted_faces.begin(),
                       ghosted_faces.end());
      Table<2, unsigned int> cell_and_face_to_faces(
        (row_starts.size() - 1) / start_components.back(),
        2 * shape_info(0, 0).n_dimensions);
      cell_and_face_to_faces.fill(numbers::invalid_unsigned_int);
      for (unsigned int f = 0; f < all_faces.size(); ++f)
        {
          cell_and_face_to_faces(all_faces[f].cells_interior[0],
                                 all_faces[f].interior_face_no) = f;
          Assert(all_faces[f].cells_exterior[0] !=
                   numbers::invalid_unsigned_int,
                 ExcInternalError());
          cell_and_face_to_faces(all_faces[f].cells_exterior[0],
                                 all_faces[f].exterior_face_no) = f;
        }

      // lambda function to detect objects on face pairs
      const auto loop_over_faces =
        [&](const std::function<
            void(const unsigned int, const unsigned int, const bool)> &fu) {
          for (const auto &face : inner_faces)
            {
              AssertIndexRange(face.cells_interior[0], n_owned_cells);
              fu(face.cells_exterior[0], face.exterior_face_no, false /*flag*/);
            }
        };

      const auto loop_over_all_faces =
        [&](const std::function<
            void(const unsigned int, const unsigned int, const bool)> &fu) {
          for (unsigned int c = 0; c < cell_and_face_to_faces.size(0); ++c)
            for (unsigned int d = 0; d < cell_and_face_to_faces.size(1); ++d)
              {
                const unsigned int f = cell_and_face_to_faces(c, d);
                if (f == numbers::invalid_unsigned_int)
                  continue;

                const unsigned int cell_m = all_faces[f].cells_interior[0];
                const unsigned int cell_p = all_faces[f].cells_exterior[0];

                const bool ext = c == cell_m;

                if (ext && cell_p == numbers::invalid_unsigned_int)
                  continue;

                const unsigned int p       = ext ? cell_p : cell_m;
                const unsigned int face_no = ext ?
                                               all_faces[f].exterior_face_no :
                                               all_faces[f].interior_face_no;

                fu(p, face_no, true);
              }
        };

      const auto process_values =
        [&](
          std::shared_ptr<const Utilities::MPI::Partitioner>
            &vector_partitioner_values,
          const std::function<void(
            const std::function<void(
              const unsigned int, const unsigned int, const bool)> &)> &loop) {
          bool all_nodal_and_tensorial = shape_info.size(1) == 1;

          if (all_nodal_and_tensorial)
            for (unsigned int c = 0; c < n_base_elements; ++c)
              {
                const auto &si =
                  shape_info(global_base_element_offset + c, 0).data.front();
                if (!si.nodal_at_cell_boundaries ||
                    (si.element_type ==
                     MatrixFreeFunctions::ElementType::tensor_none))
                  all_nodal_and_tensorial = false;
              }

          if (all_nodal_and_tensorial == false)
            vector_partitioner_values = vector_partitioner;
          else
            {
              bool has_noncontiguous_cell = false;

              loop([&](const unsigned int cell_no,
                       const unsigned int face_no,
                       const bool         flag) {
                const unsigned int index =
                  dof_indices_contiguous[dof_access_cell][cell_no];
                if (flag || (index != numbers::invalid_unsigned_int &&
                             index >= part.locally_owned_size()))
                  {
                    const unsigned int stride =
                      dof_indices_interleave_strides[dof_access_cell][cell_no];
                    unsigned int i = 0;
                    for (unsigned int e = 0; e < n_base_elements; ++e)
                      for (unsigned int c = 0; c < n_components[e]; ++c)
                        {
                          const ShapeInfo<double> &shape =
                            shape_info(global_base_element_offset + e, 0);
                          for (unsigned int j = 0;
                               j < shape.dofs_per_component_on_face;
                               ++j)
                            ghost_indices.push_back(part.local_to_global(
                              index + i +
                              shape.face_to_cell_index_nodal(face_no, j) *
                                stride));
                          i += shape.dofs_per_component_on_cell * stride;
                        }
                    AssertDimension(i, dofs_per_cell[0] * stride);
                  }
                else if (index == numbers::invalid_unsigned_int)
                  has_noncontiguous_cell = true;
              });
              has_noncontiguous_cell =
                Utilities::MPI::min<int>(static_cast<int>(
                                           has_noncontiguous_cell),
                                         part.get_mpi_communicator()) != 0;

              std::sort(ghost_indices.begin(), ghost_indices.end());
              IndexSet compressed_set(part.size());
              compressed_set.add_indices(ghost_indices.begin(),
                                         ghost_indices.end());
              compressed_set.subtract_set(part.locally_owned_range());
              const bool all_ghosts_equal =
                Utilities::MPI::min<int>(static_cast<int>(
                                           compressed_set.n_elements() ==
                                           part.ghost_indices().n_elements()),
                                         part.get_mpi_communicator()) != 0;
              if (all_ghosts_equal || has_noncontiguous_cell)
                vector_partitioner_values = vector_partitioner;
              else
                {
                  vector_partitioner_values =
                    std::make_shared<Utilities::MPI::Partitioner>(
                      part.locally_owned_range(), part.get_mpi_communicator());
                  const_cast<Utilities::MPI::Partitioner *>(
                    vector_partitioner_values.get())
                    ->set_ghost_indices(compressed_set, part.ghost_indices());
                }
            }
        };


      const auto process_gradients =
        [&](
          const std::shared_ptr<const Utilities::MPI::Partitioner>
            &vector_partitoner_values,
          std::shared_ptr<const Utilities::MPI::Partitioner>
            &vector_partitioner_gradients,
          const std::function<void(
            const std::function<void(
              const unsigned int, const unsigned int, const bool)> &)> &loop) {
          bool all_hermite = shape_info.size(1) == 1;

          if (all_hermite)
            for (unsigned int c = 0; c < n_base_elements; ++c)
              if (shape_info(global_base_element_offset + c, 0).element_type !=
                  MatrixFreeFunctions::tensor_symmetric_hermite)
                all_hermite = false;
          if (all_hermite == false ||
              vector_partitoner_values.get() == vector_partitioner.get())
            vector_partitioner_gradients = vector_partitioner;
          else
            {
              loop([&](const unsigned int cell_no,
                       const unsigned int face_no,
                       const bool         flag) {
                const unsigned int index =
                  dof_indices_contiguous[dof_access_cell][cell_no];
                if (flag || (index != numbers::invalid_unsigned_int &&
                             index >= part.locally_owned_size()))
                  {
                    const unsigned int stride =
                      dof_indices_interleave_strides[dof_access_cell][cell_no];
                    unsigned int i = 0;
                    for (unsigned int e = 0; e < n_base_elements; ++e)
                      for (unsigned int c = 0; c < n_components[e]; ++c)
                        {
                          const ShapeInfo<double> &shape =
                            shape_info(global_base_element_offset + e, 0);
                          for (unsigned int j = 0;
                               j < 2 * shape.dofs_per_component_on_face;
                               ++j)
                            ghost_indices.push_back(part.local_to_global(
                              index + i +
                              shape.face_to_cell_index_hermite(face_no, j) *
                                stride));
                          i += shape.dofs_per_component_on_cell * stride;
                        }
                    AssertDimension(i, dofs_per_cell[0] * stride);
                  }
              });
              std::sort(ghost_indices.begin(), ghost_indices.end());
              IndexSet compressed_set(part.size());
              compressed_set.add_indices(ghost_indices.begin(),
                                         ghost_indices.end());
              compressed_set.subtract_set(part.locally_owned_range());
              const bool all_ghosts_equal =
                Utilities::MPI::min<int>(static_cast<int>(
                                           compressed_set.n_elements() ==
                                           part.ghost_indices().n_elements()),
                                         part.get_mpi_communicator()) != 0;
              if (all_ghosts_equal)
                vector_partitioner_gradients = vector_partitioner;
              else
                {
                  vector_partitioner_gradients =
                    std::make_shared<Utilities::MPI::Partitioner>(
                      part.locally_owned_range(), part.get_mpi_communicator());
                  const_cast<Utilities::MPI::Partitioner *>(
                    vector_partitioner_gradients.get())
                    ->set_ghost_indices(compressed_set, part.ghost_indices());
                }
            }
        };

      std::shared_ptr<const Utilities::MPI::Partitioner> temp_1, temp_2, temp_3,
        temp_4;

      // partitioner 1: values on faces
      process_values(temp_1, loop_over_faces);

      // partitioner 2: values and gradients on faces
      process_gradients(temp_1, temp_2, loop_over_faces);

      if (fill_cell_centric)
        {
          ghost_indices.clear();
          // partitioner 3: values on all faces
          process_values(temp_3, loop_over_all_faces);
          // partitioner 4: values and gradients on faces
          process_gradients(temp_3, temp_4, loop_over_all_faces);
        }
      else
        {
          temp_3 = std::make_shared<Utilities::MPI::Partitioner>(
            part.locally_owned_range(), part.get_mpi_communicator());
          temp_4 = std::make_shared<Utilities::MPI::Partitioner>(
            part.locally_owned_range(), part.get_mpi_communicator());
        }

      if (use_vector_data_exchanger_full == false)
        {
          vector_exchanger_face_variants[1] = std::make_shared<
            MatrixFreeFunctions::VectorDataExchange::PartitionerWrapper>(
            temp_1);
          vector_exchanger_face_variants[2] = std::make_shared<
            MatrixFreeFunctions::VectorDataExchange::PartitionerWrapper>(
            temp_2);
          vector_exchanger_face_variants[3] = std::make_shared<
            MatrixFreeFunctions::VectorDataExchange::PartitionerWrapper>(
            temp_3);
          vector_exchanger_face_variants[4] = std::make_shared<
            MatrixFreeFunctions::VectorDataExchange::PartitionerWrapper>(
            temp_4);
        }
      else
        {
          vector_exchanger_face_variants[1] =
            std::make_shared<MatrixFreeFunctions::VectorDataExchange::Full>(
              temp_1, communicator_sm);
          vector_exchanger_face_variants[2] =
            std::make_shared<MatrixFreeFunctions::VectorDataExchange::Full>(
              temp_2, communicator_sm);
          vector_exchanger_face_variants[3] =
            std::make_shared<MatrixFreeFunctions::VectorDataExchange::Full>(
              temp_3, communicator_sm);
          vector_exchanger_face_variants[4] =
            std::make_shared<MatrixFreeFunctions::VectorDataExchange::Full>(
              temp_4, communicator_sm);
        }
    }



    void
    DoFInfo::compute_shared_memory_contiguous_indices(
      std::array<std::vector<std::pair<unsigned int, unsigned int>>, 3>
        &cell_indices_contiguous_sm)
    {
      AssertDimension(dofs_per_cell.size(), 1);

      for (unsigned int i = 0; i < 3; ++i)
        {
          dof_indices_contiguous_sm[i].resize(
            cell_indices_contiguous_sm[i].size());

          for (unsigned int j = 0; j < cell_indices_contiguous_sm[i].size();
               ++j)
            if (cell_indices_contiguous_sm[i][j].first !=
                numbers::invalid_unsigned_int)
              dof_indices_contiguous_sm[i][j] = {
                cell_indices_contiguous_sm[i][j].first,
                cell_indices_contiguous_sm[i][j].second * dofs_per_cell[0]};
            else
              dof_indices_contiguous_sm[i][j] = {numbers::invalid_unsigned_int,
                                                 numbers::invalid_unsigned_int};
        }
    }



    namespace internal
    {
      // We construct the connectivity graph in parallel. we use one lock for
      // 256 degrees of freedom to keep the number of locks down to a
      // reasonable level and reduce the cost of locking to some extent.
      static constexpr unsigned int bucket_size_threading = 256;



      void
      compute_row_lengths(const unsigned int         begin,
                          const unsigned int         end,
                          const DoFInfo             &dof_info,
                          std::vector<std::mutex>   &mutexes,
                          std::vector<unsigned int> &row_lengths)
      {
        std::vector<unsigned int> scratch;
        const unsigned int n_components = dof_info.start_components.back();
        for (unsigned int block = begin; block < end; ++block)
          {
            scratch.assign(
              dof_info.dof_indices.data() +
                dof_info.row_starts[block * n_components].first,
              dof_info.dof_indices.data() +
                dof_info.row_starts[(block + 1) * n_components].first);
            std::sort(scratch.begin(), scratch.end());

            const std::vector<unsigned int>::const_iterator end_unique =
              std::unique(scratch.begin(), scratch.end());
            for (std::vector<unsigned int>::const_iterator it = scratch.begin();
                 it != end_unique;
                 /* update in loop body */)
              {
                // In this code, the procedure is that we insert all elements
                // that are within the range of one lock at once
                const unsigned int next_bucket =
                  (*it / bucket_size_threading + 1) * bucket_size_threading;

                std::lock_guard<std::mutex> lock(
                  mutexes[*it / bucket_size_threading]);
                for (; it != end_unique && *it < next_bucket; ++it)
                  {
                    AssertIndexRange(*it, row_lengths.size());
                    ++row_lengths[*it];
                  }
              }
          }
      }

      void
      fill_connectivity_dofs(const unsigned int               begin,
                             const unsigned int               end,
                             const DoFInfo                   &dof_info,
                             const std::vector<unsigned int> &row_lengths,
                             std::vector<std::mutex>         &mutexes,
                             dealii::SparsityPattern         &connectivity_dof)
      {
        std::vector<unsigned int> scratch;
        const unsigned int n_components = dof_info.start_components.back();
        for (unsigned int block = begin; block < end; ++block)
          {
            scratch.assign(
              dof_info.dof_indices.data() +
                dof_info.row_starts[block * n_components].first,
              dof_info.dof_indices.data() +
                dof_info.row_starts[(block + 1) * n_components].first);
            std::sort(scratch.begin(), scratch.end());

            const std::vector<unsigned int>::const_iterator end_unique =
              std::unique(scratch.begin(), scratch.end());
            for (std::vector<unsigned int>::const_iterator it = scratch.begin();
                 it != end_unique;
                 /* update in loop body */)
              {
                const unsigned int next_bucket =
                  (*it / bucket_size_threading + 1) * bucket_size_threading;

                std::lock_guard<std::mutex> lock(
                  mutexes[*it / bucket_size_threading]);
                for (; it != end_unique && *it < next_bucket; ++it)
                  if (row_lengths[*it] > 0)
                    connectivity_dof.add(*it, block);
              }
          }
      }



      void
      fill_connectivity(const unsigned int               begin,
                        const unsigned int               end,
                        const DoFInfo                   &dof_info,
                        const std::vector<unsigned int> &renumbering,
                        const dealii::SparsityPattern   &connectivity_dof,
                        DynamicSparsityPattern          &connectivity)
      {
        ordered_vector     row_entries;
        const unsigned int n_components = dof_info.start_components.back();
        for (unsigned int block = begin; block < end; ++block)
          {
            row_entries.clear();

            const unsigned int
              *it = dof_info.dof_indices.data() +
                    dof_info.row_starts[block * n_components].first,
              *end_cell = dof_info.dof_indices.data() +
                          dof_info.row_starts[(block + 1) * n_components].first;
            for (; it != end_cell; ++it)
              {
                SparsityPattern::iterator sp = connectivity_dof.begin(*it);
                std::vector<types::global_dof_index>::iterator insert_pos =
                  row_entries.begin();
                for (; sp != connectivity_dof.end(*it); ++sp)
                  if (sp->column() != block)
                    row_entries.insert(renumbering[sp->column()], insert_pos);
              }
            connectivity.add_entries(renumbering[block],
                                     row_entries.begin(),
                                     row_entries.end());
          }
      }

    } // namespace internal

    void
    DoFInfo::make_connectivity_graph(
      const TaskInfo                  &task_info,
      const std::vector<unsigned int> &renumbering,
      DynamicSparsityPattern          &connectivity) const
    {
      unsigned int n_rows = (vector_partitioner->local_range().second -
                             vector_partitioner->local_range().first) +
                            vector_partitioner->ghost_indices().n_elements();

      // Avoid square sparsity patterns that allocate the diagonal entry
      if (n_rows == task_info.n_active_cells)
        ++n_rows;

      // first determine row lengths
      std::vector<unsigned int> row_lengths(n_rows);
      std::vector<std::mutex> mutexes(n_rows / internal::bucket_size_threading +
                                      1);
      dealii::parallel::apply_to_subranges(
        0,
        task_info.n_active_cells,
        [this, &mutexes, &row_lengths](const unsigned int begin,
                                       const unsigned int end) {
          internal::compute_row_lengths(
            begin, end, *this, mutexes, row_lengths);
        },
        20);

      // disregard dofs that only sit on a single cell because they cannot
      // couple
      for (unsigned int row = 0; row < n_rows; ++row)
        if (row_lengths[row] <= 1)
          row_lengths[row] = 0;

      // Create a temporary sparsity pattern that holds to each degree of
      // freedom on which cells it appears, i.e., store the connectivity
      // between cells and dofs
      SparsityPattern connectivity_dof(n_rows,
                                       task_info.n_active_cells,
                                       row_lengths);
      dealii::parallel::apply_to_subranges(
        0,
        task_info.n_active_cells,
        [this, &row_lengths, &mutexes, &connectivity_dof](
          const unsigned int begin, const unsigned int end) {
          internal::fill_connectivity_dofs(
            begin, end, *this, row_lengths, mutexes, connectivity_dof);
        },
        20);
      connectivity_dof.compress();


      // Invert renumbering for use in fill_connectivity.
      std::vector<unsigned int> reverse_numbering(task_info.n_active_cells);
      reverse_numbering = Utilities::invert_permutation(renumbering);

      // From the above connectivity between dofs and cells, we can finally
      // create a connectivity list between cells. The connectivity graph
      // should apply the renumbering, i.e., the entry for cell j is the entry
      // for cell renumbering[j] in the original ordering.
      dealii::parallel::apply_to_subranges(
        0,
        task_info.n_active_cells,
        [this, &reverse_numbering, &connectivity_dof, &connectivity](
          const unsigned int begin, const unsigned int end) {
          internal::fill_connectivity(begin,
                                      end,
                                      *this,
                                      reverse_numbering,
                                      connectivity_dof,
                                      connectivity);
        },
        20);
    }



    void
    DoFInfo::compute_dof_renumbering(
      std::vector<types::global_dof_index> &renumbering)
    {
      const unsigned int locally_owned_size =
        vector_partitioner->locally_owned_size();
      renumbering.resize(0);
      renumbering.resize(locally_owned_size, numbers::invalid_dof_index);

      types::global_dof_index counter      = 0;
      const unsigned int      n_components = start_components.back();
      const unsigned int      n_cell_batches =
        n_vectorization_lanes_filled[dof_access_cell].size();
      Assert(n_cell_batches <=
               (row_starts.size() - 1) / vectorization_length / n_components,
             ExcInternalError());
      for (unsigned int cell_no = 0; cell_no < n_cell_batches; ++cell_no)
        {
          // do not renumber in case we have constraints
          if (row_starts[cell_no * n_components * vectorization_length]
                .second ==
              row_starts[(cell_no + 1) * n_components * vectorization_length]
                .second)
            {
              const unsigned int ndofs =
                dofs_per_cell.size() == 1 ?
                  dofs_per_cell[0] :
                  (dofs_per_cell[cell_active_fe_index.size() > 0 ?
                                   cell_active_fe_index[cell_no] :
                                   0]);
              const unsigned int *dof_ind =
                dof_indices.data() +
                row_starts[cell_no * n_components * vectorization_length].first;
              for (unsigned int i = 0; i < ndofs; ++i)
                for (unsigned int j = 0;
                     j < n_vectorization_lanes_filled[dof_access_cell][cell_no];
                     ++j)
                  if (dof_ind[j * ndofs + i] < locally_owned_size)
                    if (renumbering[dof_ind[j * ndofs + i]] ==
                        numbers::invalid_dof_index)
                      renumbering[dof_ind[j * ndofs + i]] = counter++;
            }
        }

      AssertIndexRange(counter, locally_owned_size + 1);
      for (types::global_dof_index &dof_index : renumbering)
        if (dof_index == numbers::invalid_dof_index)
          dof_index = counter++;

      // transform indices to global index space
      for (types::global_dof_index &dof_index : renumbering)
        dof_index = vector_partitioner->local_to_global(dof_index);

      AssertDimension(counter, renumbering.size());
    }



    std::size_t
    DoFInfo::memory_consumption() const
    {
      std::size_t memory = sizeof(*this);
      for (const auto &storage : index_storage_variants)
        memory += storage.capacity() * sizeof(storage[0]);
      memory +=
        (row_starts.capacity() * sizeof(std::pair<unsigned int, unsigned int>));
      memory += MemoryConsumption::memory_consumption(dof_indices);
      memory += MemoryConsumption::memory_consumption(dof_indices_interleaved);
      memory += MemoryConsumption::memory_consumption(dof_indices_contiguous);
      memory +=
        MemoryConsumption::memory_consumption(dof_indices_contiguous_sm);
      memory +=
        MemoryConsumption::memory_consumption(dof_indices_interleave_strides);
      memory +=
        MemoryConsumption::memory_consumption(n_vectorization_lanes_filled);
      memory += MemoryConsumption::memory_consumption(
        hanging_node_constraint_masks_comp);
      memory +=
        MemoryConsumption::memory_consumption(hanging_node_constraint_masks);
      memory += MemoryConsumption::memory_consumption(constrained_dofs);
      memory += MemoryConsumption::memory_consumption(row_starts_plain_indices);
      memory += MemoryConsumption::memory_consumption(plain_dof_indices);
      memory += MemoryConsumption::memory_consumption(constraint_indicator);
      memory += MemoryConsumption::memory_consumption(*vector_partitioner);
      memory += MemoryConsumption::memory_consumption(n_components);
      memory += MemoryConsumption::memory_consumption(start_components);
      memory += MemoryConsumption::memory_consumption(component_to_base_index);
      memory +=
        MemoryConsumption::memory_consumption(component_dof_indices_offset);
      memory += MemoryConsumption::memory_consumption(dofs_per_cell);
      memory += MemoryConsumption::memory_consumption(dofs_per_face);
      memory += MemoryConsumption::memory_consumption(cell_active_fe_index);
      memory += MemoryConsumption::memory_consumption(fe_index_conversion);
      memory +=
        MemoryConsumption::memory_consumption(vector_zero_range_list_index);
      memory += MemoryConsumption::memory_consumption(vector_zero_range_list);
      memory += MemoryConsumption::memory_consumption(cell_loop_pre_list_index);
      memory += MemoryConsumption::memory_consumption(cell_loop_pre_list);
      memory +=
        MemoryConsumption::memory_consumption(cell_loop_post_list_index);
      memory += MemoryConsumption::memory_consumption(cell_loop_post_list);
      return memory;
    }
  } // namespace MatrixFreeFunctions
} // namespace internal

namespace internal
{
  namespace MatrixFreeFunctions
  {
    template void
    DoFInfo::read_dof_indices<double>(
      const std::vector<types::global_dof_index> &,
      const std::vector<types::global_dof_index> &,
      const bool,
      const dealii::AffineConstraints<double> &,
      const unsigned int,
      ConstraintValues<double> &,
      bool &);

    template void
    DoFInfo::read_dof_indices<float>(
      const std::vector<types::global_dof_index> &,
      const std::vector<types::global_dof_index> &,
      const bool,
      const dealii::AffineConstraints<float> &,
      const unsigned int,
      ConstraintValues<double> &,
      bool &);

    template bool
    DoFInfo::process_hanging_node_constraints<1>(
      const HangingNodes<1> &,
      const std::vector<std::vector<unsigned int>> &,
      const unsigned int,
      const TriaIterator<DoFCellAccessor<1, 1, false>> &,
      std::vector<types::global_dof_index> &);
    template bool
    DoFInfo::process_hanging_node_constraints<2>(
      const HangingNodes<2> &,
      const std::vector<std::vector<unsigned int>> &,
      const unsigned int,
      const TriaIterator<DoFCellAccessor<2, 2, false>> &,
      std::vector<types::global_dof_index> &);
    template bool
    DoFInfo::process_hanging_node_constraints<3>(
      const HangingNodes<3> &,
      const std::vector<std::vector<unsigned int>> &,
      const unsigned int,
      const TriaIterator<DoFCellAccessor<3, 3, false>> &,
      std::vector<types::global_dof_index> &);

    template void
    DoFInfo::compute_face_index_compression<1>(
      const std::vector<FaceToCellTopology<1>> &,
      bool);
    template void
    DoFInfo::compute_face_index_compression<2>(
      const std::vector<FaceToCellTopology<2>> &,
      bool);
    template void
    DoFInfo::compute_face_index_compression<4>(
      const std::vector<FaceToCellTopology<4>> &,
      bool);
    template void
    DoFInfo::compute_face_index_compression<8>(
      const std::vector<FaceToCellTopology<8>> &,
      bool);
    template void
    DoFInfo::compute_face_index_compression<16>(
      const std::vector<FaceToCellTopology<16>> &,
      bool);

    template void
    DoFInfo::compute_vector_zero_access_pattern<1>(
      const TaskInfo &,
      const std::vector<FaceToCellTopology<1>> &);
    template void
    DoFInfo::compute_vector_zero_access_pattern<2>(
      const TaskInfo &,
      const std::vector<FaceToCellTopology<2>> &);
    template void
    DoFInfo::compute_vector_zero_access_pattern<4>(
      const TaskInfo &,
      const std::vector<FaceToCellTopology<4>> &);
    template void
    DoFInfo::compute_vector_zero_access_pattern<8>(
      const TaskInfo &,
      const std::vector<FaceToCellTopology<8>> &);
    template void
    DoFInfo::compute_vector_zero_access_pattern<16>(
      const TaskInfo &,
      const std::vector<FaceToCellTopology<16>> &);

    template void
    DoFInfo::print_memory_consumption<std::ostream>(std::ostream &,
                                                    const TaskInfo &) const;
    template void
    DoFInfo::print_memory_consumption<ConditionalOStream>(
      ConditionalOStream &,
      const TaskInfo &) const;

    template void
    DoFInfo::print<double>(const std::vector<double> &,
                           const std::vector<unsigned int> &,
                           std::ostream &) const;

    template void
    DoFInfo::print<float>(const std::vector<float> &,
                          const std::vector<unsigned int> &,
                          std::ostream &) const;
  } // namespace MatrixFreeFunctions
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
