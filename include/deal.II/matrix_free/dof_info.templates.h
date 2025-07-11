// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_matrix_free_dof_info_templates_h
#define dealii_matrix_free_dof_info_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/floating_point_comparator.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/parallel.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/constraint_info.h>
#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/hanging_nodes_internal.h>
#include <deal.II/matrix_free/task_info.h>

#include <limits>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    template <typename number>
    void
    DoFInfo::read_dof_indices(
      const std::vector<types::global_dof_index> &local_indices_resolved,
      const std::vector<types::global_dof_index> &local_indices,
      const bool                                  cell_has_hanging_nodes,
      const dealii::AffineConstraints<number>    &constraints,
      const unsigned int                          cell_number,
      ConstraintValues<double>                   &constraint_values,
      bool                                       &cell_at_subdomain_boundary)
    {
      Assert(vector_partitioner.get() != nullptr, ExcInternalError());
      const types::global_dof_index first_owned =
        vector_partitioner->local_range().first;
      const types::global_dof_index last_owned =
        vector_partitioner->local_range().second;
      Assert(last_owned - first_owned <
               std::numeric_limits<unsigned int>::max(),
             ExcMessage("The size local range of owned indices must not "
                        "exceed the size of unsigned int"));
      const unsigned int n_owned = last_owned - first_owned;

      Assert(dofs_per_cell.size() == 1 ||
               cell_number < cell_active_fe_index.size(),
             ExcInternalError());
      const unsigned int fe_index =
        dofs_per_cell.size() == 1 ? 0 : cell_active_fe_index[cell_number];
      const unsigned int dofs_this_cell = dofs_per_cell[fe_index];
      const unsigned int n_components   = start_components.back();
      unsigned int       i              = 0;
      for (unsigned int comp = 0; comp < n_components; ++comp)
        {
          std::pair<unsigned short, unsigned short> constraint_iterator(0, 0);
          const auto next = component_dof_indices_offset[fe_index][comp + 1];
          for (; i < next; ++i)
            {
              types::global_dof_index current_dof = local_indices_resolved[i];
              const auto             *entries_ptr =
                constraints.get_constraint_entries(current_dof);

              // dof is constrained
              if (entries_ptr != nullptr)
                {
                  // in case we want to access plain indices, we need to know
                  // about the location of constrained indices as well (all the
                  // other indices are collected by the cases below)
                  if (current_dof < first_owned || current_dof >= last_owned)
                    {
                      ghost_dofs.push_back(current_dof);
                      cell_at_subdomain_boundary = true;
                    }

                  // check whether this dof is identity constrained to another
                  // dof. then we can simply insert that dof and there is no
                  // need to actually resolve the constraint entries
                  const auto                   &entries   = *entries_ptr;
                  const types::global_dof_index n_entries = entries.size();
                  if (n_entries == 1 &&
                      std::abs(entries[0].second - 1.) <
                        100 * std::numeric_limits<double>::epsilon())
                    {
                      current_dof = entries[0].first;
                      goto no_constraint;
                    }

                  // append a new index to the indicators
                  constraint_indicator.push_back(constraint_iterator);
                  constraint_indicator.back().second =
                    constraint_values.insert_entries(entries);

                  // reset constraint iterator for next round
                  constraint_iterator.first = 0;

                  // add the local_to_global indices computed in the
                  // insert_entries function. transform the index to local index
                  // space or mark it as ghost if necessary
                  if (n_entries > 0)
                    {
                      const std::vector<types::global_dof_index>
                        &constraint_indices =
                          constraint_values.constraint_indices;
                      for (unsigned int j = 0; j < n_entries; ++j)
                        {
                          if (constraint_indices[j] < first_owned ||
                              constraint_indices[j] >= last_owned)
                            {
                              dof_indices.push_back(n_owned +
                                                    ghost_dofs.size());

                              // collect ghosts so that we can later construct
                              // an IndexSet for them. also store whether the
                              // current cell is on the boundary
                              ghost_dofs.push_back(constraint_indices[j]);
                              cell_at_subdomain_boundary = true;
                            }
                          else
                            // not ghost, so transform to the local index space
                            // directly
                            dof_indices.push_back(static_cast<unsigned int>(
                              constraint_indices[j] - first_owned));
                        }
                    }
                }
              else
                {
                no_constraint:
                  // Not constrained, we simply have to add the local index to
                  // the indices_local_to_global list and increment constraint
                  // iterator. transform to local index space/mark as ghost
                  if (current_dof < first_owned || current_dof >= last_owned)
                    {
                      ghost_dofs.push_back(current_dof);
                      current_dof = n_owned + ghost_dofs.size() - 1;
                      cell_at_subdomain_boundary = true;
                    }
                  else
                    current_dof -= first_owned;

                  dof_indices.push_back(static_cast<unsigned int>(current_dof));

                  // make sure constraint_iterator.first is always within the
                  // bounds of unsigned short
                  Assert(constraint_iterator.first <
                           (1 << (8 * sizeof(unsigned short))) - 1,
                         ExcInternalError());
                  constraint_iterator.first++;
                }
            }
          row_starts[cell_number * n_components + comp + 1].first =
            dof_indices.size();
          row_starts[cell_number * n_components + comp + 1].second =
            constraint_indicator.size();
        }

      // now to the plain indices: in case we have constraints on this cell,
      // store the indices without the constraints resolve once again
      if (store_plain_indices == true)
        {
          if (cell_number == 0)
            row_starts_plain_indices.resize(
              (row_starts.size() - 1) / n_components + 1);
          row_starts_plain_indices[cell_number] = plain_dof_indices.size();
          const bool cell_has_constraints =
            cell_has_hanging_nodes ||
            (row_starts[(cell_number + 1) * n_components].second >
             row_starts[cell_number * n_components].second);
          if (cell_has_constraints == true)
            {
              for (unsigned int i = 0; i < dofs_this_cell; ++i)
                {
                  types::global_dof_index current_dof = local_indices[i];
                  if (current_dof < first_owned || current_dof >= last_owned)
                    {
                      ghost_dofs.push_back(current_dof);
                      current_dof = n_owned + ghost_dofs.size() - 1;
                      cell_at_subdomain_boundary = true;
                    }
                  else
                    current_dof -= first_owned;
                  plain_dof_indices.push_back(
                    static_cast<unsigned int>(current_dof));
                }
            }
        }
    }



    template <int dim>
    bool
    DoFInfo::process_hanging_node_constraints(
      const HangingNodes<dim>                      &hanging_nodes,
      const std::vector<std::vector<unsigned int>> &lexicographic_mapping,
      const unsigned int                            cell_number,
      const TriaIterator<DoFCellAccessor<dim, dim, false>> &cell,
      std::vector<types::global_dof_index>                 &dof_indices)
    {
      if (this->hanging_node_constraint_masks_comp.empty())
        return false;

      // 2) determine the refinement configuration of the cell
      const auto refinement_configuration =
        hanging_nodes.compute_refinement_configuration(cell);

      if (refinement_configuration == ConstraintKinds::unconstrained)
        return false;

      // 3) update DoF indices of cell for specified components
      hanging_nodes.update_dof_indices(cell,
                                       {},
                                       lexicographic_mapping,
                                       hanging_node_constraint_masks_comp,
                                       refinement_configuration,
                                       dof_indices);

      hanging_node_constraint_masks[cell_number] =
        compress(refinement_configuration, dim);

      return true;
    }



    template <int length>
    void
    DoFInfo::compute_face_index_compression(
      const std::vector<FaceToCellTopology<length>> &faces,
      const bool hold_all_faces_to_owned_cells)
    {
      AssertDimension(length, vectorization_length);

      index_storage_variants[dof_access_face_interior].resize(
        faces.size(), IndexStorageVariants::full);
      dof_indices_contiguous[dof_access_face_interior].resize(
        faces.size() * length, numbers::invalid_unsigned_int);
      dof_indices_interleave_strides[dof_access_face_interior].resize(
        faces.size() * length, numbers::invalid_unsigned_int);
      n_vectorization_lanes_filled[dof_access_face_interior].resize(
        faces.size());

      // all inner faces come before the boundary faces
      unsigned int n_inner_faces = 0;
      for (; n_inner_faces < faces.size(); ++n_inner_faces)
        if (faces[n_inner_faces].cells_exterior[0] ==
            numbers::invalid_unsigned_int)
          break;

      // all boundary faces come after the inner faces and before the ghosted
      // inner faces
      unsigned int n_boundary_faces = 0;
      for (; n_inner_faces + n_boundary_faces < faces.size();
           ++n_boundary_faces)
        if (faces[n_inner_faces + n_boundary_faces].cells_exterior[0] !=
            numbers::invalid_unsigned_int)
          break;

      const unsigned int size_exterior_faces =
        hold_all_faces_to_owned_cells ? faces.size() : n_inner_faces;

      index_storage_variants[dof_access_face_exterior].resize(
        size_exterior_faces, IndexStorageVariants::full);
      dof_indices_contiguous[dof_access_face_exterior].resize(
        size_exterior_faces * length, numbers::invalid_unsigned_int);
      dof_indices_interleave_strides[dof_access_face_exterior].resize(
        faces.size() * length, numbers::invalid_unsigned_int);
      n_vectorization_lanes_filled[dof_access_face_exterior].resize(
        size_exterior_faces);

      for (unsigned int face = 0; face < faces.size(); ++face)
        {
          auto face_computation = [&](const DoFAccessIndex face_index,
                                      const std::array<unsigned int, length>
                                        &cell_indices_face) {
            bool is_contiguous      = false;
            bool is_interleaved     = false;
            bool needs_full_storage = false;
            for (unsigned int v = 0;
                 v < length &&
                 cell_indices_face[v] != numbers::invalid_unsigned_int;
                 ++v)
              {
                n_vectorization_lanes_filled[face_index][face]++;
                if (index_storage_variants[dof_access_cell]
                                          [cell_indices_face[v] / length] >=
                    IndexStorageVariants::interleaved_contiguous)
                  is_interleaved = true;
                if (index_storage_variants[dof_access_cell]
                                          [cell_indices_face[v] / length] ==
                    IndexStorageVariants::contiguous)
                  is_contiguous = true;
                if (index_storage_variants[dof_access_cell]
                                          [cell_indices_face[v] / length] >=
                    IndexStorageVariants::contiguous)
                  dof_indices_interleave_strides[face_index][face * length +
                                                             v] =
                    dof_indices_interleave_strides[dof_access_cell]
                                                  [cell_indices_face[v]];
                if (index_storage_variants[dof_access_cell]
                                          [cell_indices_face[v] / length] <
                    IndexStorageVariants::contiguous)
                  needs_full_storage = true;
              }
            Assert(!(is_interleaved && is_contiguous),
                   ExcMessage("Unsupported index compression found"));

            if (is_interleaved || is_contiguous)
              for (unsigned int v = 0;
                   v < n_vectorization_lanes_filled[face_index][face];
                   ++v)
                dof_indices_contiguous[face_index][face * length + v] =
                  dof_indices_contiguous[dof_access_cell][cell_indices_face[v]];
            if (is_interleaved)
              {
                bool is_also_contiguous =
                  n_vectorization_lanes_filled[face_index][face] == length;
                for (unsigned int v = 0;
                     v < n_vectorization_lanes_filled[face_index][face];
                     ++v)
                  if (dof_indices_contiguous[face_index][face * length + v] !=
                        dof_indices_contiguous[face_index][face * length] + v ||
                      dof_indices_interleave_strides[dof_access_cell]
                                                    [cell_indices_face[v]] !=
                        length)
                    is_also_contiguous = false;

                if (is_also_contiguous)
                  {
                    index_storage_variants[face_index][face] =
                      IndexStorageVariants::interleaved_contiguous;
                  }
                else
                  {
                    bool all_indices_same_offset =
                      n_vectorization_lanes_filled[face_index][face] == length;
                    for (unsigned int v = 0;
                         v < n_vectorization_lanes_filled[face_index][face];
                         ++v)
                      if (dof_indices_interleave_strides
                            [dof_access_cell][cell_indices_face[v]] != length)
                        all_indices_same_offset = false;
                    if (all_indices_same_offset)
                      index_storage_variants[face_index][face] =
                        IndexStorageVariants::interleaved_contiguous_strided;
                    else
                      index_storage_variants[face_index][face] =
                        IndexStorageVariants::
                          interleaved_contiguous_mixed_strides;
                  }
              }
            else if (is_contiguous && !needs_full_storage)
              index_storage_variants[face_index][face] =
                IndexStorageVariants::contiguous;
            else
              index_storage_variants[face_index][face] =
                IndexStorageVariants::full;
          };

          face_computation(dof_access_face_interior,
                           faces[face].cells_interior);
          if (face < n_inner_faces ||
              (hold_all_faces_to_owned_cells &&
               face >= (n_inner_faces + n_boundary_faces)))
            face_computation(dof_access_face_exterior,
                             faces[face].cells_exterior);
        }
    }



    template <int length>
    void
    DoFInfo::compute_vector_zero_access_pattern(
      const TaskInfo                                &task_info,
      const std::vector<FaceToCellTopology<length>> &faces)
    {
      // compute a list that tells us the first time a degree of freedom is
      // touched by a cell
      AssertDimension(length, vectorization_length);
      const unsigned int n_components = start_components.back();
      const unsigned int n_dofs = vector_partitioner->locally_owned_size() +
                                  vector_partitioner->n_ghost_indices();
      std::vector<unsigned int> touched_first_by(
        (n_dofs + chunk_size_zero_vector - 1) / chunk_size_zero_vector,
        numbers::invalid_unsigned_int);
      std::vector<unsigned int> touched_last_by(
        (n_dofs + chunk_size_zero_vector - 1) / chunk_size_zero_vector,
        numbers::invalid_unsigned_int);
      std::vector<unsigned int> cells_in_interval;
      for (unsigned int part = 0;
           part < task_info.partition_row_index.size() - 2;
           ++part)
        for (unsigned int chunk = task_info.partition_row_index[part];
             chunk < task_info.partition_row_index[part + 1];
             ++chunk)
          {
            cells_in_interval.clear();
            for (unsigned int cell = task_info.cell_partition_data[chunk];
                 cell < task_info.cell_partition_data[chunk + 1];
                 ++cell)
              for (unsigned int v = 0; v < vectorization_length; ++v)
                cells_in_interval.push_back(cell * vectorization_length + v);
            if (faces.size() > 0)
              {
                for (unsigned int face = task_info.face_partition_data[chunk];
                     face < task_info.face_partition_data[chunk + 1];
                     ++face)
                  for (unsigned int v = 0; v < vectorization_length; ++v)
                    {
                      if (faces[face].cells_interior[v] !=
                          numbers::invalid_unsigned_int)
                        cells_in_interval.push_back(
                          faces[face].cells_interior[v]);
                      if (faces[face].cells_exterior[v] !=
                          numbers::invalid_unsigned_int)
                        cells_in_interval.push_back(
                          faces[face].cells_exterior[v]);
                    }
                for (unsigned int face =
                       task_info.boundary_partition_data[chunk];
                     face < task_info.boundary_partition_data[chunk + 1];
                     ++face)
                  for (unsigned int v = 0; v < vectorization_length; ++v)
                    if (faces[face].cells_interior[v] !=
                        numbers::invalid_unsigned_int)
                      cells_in_interval.push_back(
                        faces[face].cells_interior[v]);
              }
            std::sort(cells_in_interval.begin(), cells_in_interval.end());
            cells_in_interval.erase(std::unique(cells_in_interval.begin(),
                                                cells_in_interval.end()),
                                    cells_in_interval.end());

            for (const unsigned int cell : cells_in_interval)
              {
                for (unsigned int it = row_starts[cell * n_components].first;
                     it != row_starts[(cell + 1) * n_components].first;
                     ++it)
                  {
                    const unsigned int myindex =
                      dof_indices[it] / chunk_size_zero_vector;
                    if (touched_first_by[myindex] ==
                        numbers::invalid_unsigned_int)
                      touched_first_by[myindex] = chunk;
                    touched_last_by[myindex] = chunk;
                  }
              }
          }

      // ensure that all indices are touched at least during the last round
      for (auto &index : touched_first_by)
        if (index == numbers::invalid_unsigned_int)
          index =
            task_info
              .partition_row_index[task_info.partition_row_index.size() - 2] -
            1;

      // lambda to convert from a map, with keys associated to the buckets by
      // which we sliced the index space, length chunk_size_zero_vector, and
      // values equal to the slice index which are touched by the respective
      // partition, to a "vectors-of-vectors" like data structure. Rather than
      // using the vectors, we set up a sparsity-pattern like structure where
      // one index specifies the start index (range_list_index), and the other
      // the actual ranges (range_list).
      auto convert_map_to_range_list =
        [=](const unsigned int n_partitions,
            const std::map<unsigned int, std::vector<unsigned int>> &ranges_in,
            std::vector<unsigned int> &range_list_index,
            std::vector<std::pair<unsigned int, unsigned int>> &range_list,
            const unsigned int                                  max_size) {
          range_list_index.resize(n_partitions + 1);
          range_list_index[0] = 0;
          range_list.clear();
          for (unsigned int partition = 0; partition < n_partitions;
               ++partition)
            {
              auto it = ranges_in.find(partition);
              if (it != ranges_in.end())
                {
                  for (unsigned int i = 0; i < it->second.size(); ++i)
                    {
                      const unsigned int first_i = i;
                      while (i + 1 < it->second.size() &&
                             it->second[i + 1] == it->second[i] + 1)
                        ++i;
                      range_list.emplace_back(
                        std::min(it->second[first_i] * chunk_size_zero_vector,
                                 max_size),
                        std::min((it->second[i] + 1) * chunk_size_zero_vector,
                                 max_size));
                    }
                  range_list_index[partition + 1] = range_list.size();
                }
              else
                range_list_index[partition + 1] = range_list_index[partition];
            }
        };

      // first we determine the ranges to zero the vector
      std::map<unsigned int, std::vector<unsigned int>> chunk_must_zero_vector;
      for (unsigned int i = 0; i < touched_first_by.size(); ++i)
        chunk_must_zero_vector[touched_first_by[i]].push_back(i);
      const unsigned int n_partitions =
        task_info.partition_row_index[task_info.partition_row_index.size() - 2];
      convert_map_to_range_list(n_partitions,
                                chunk_must_zero_vector,
                                vector_zero_range_list_index,
                                vector_zero_range_list,
                                vector_partitioner->locally_owned_size());

      // the other two operations only work on the local range (without
      // ghosts), so we skip the latter parts of the vector now
      touched_first_by.resize((vector_partitioner->locally_owned_size() +
                               chunk_size_zero_vector - 1) /
                              chunk_size_zero_vector);

      // set the import indices in the vector partitioner to one index higher
      // to indicate that we want to process it first. This additional index
      // is reflected in the argument 'n_partitions+1' in the
      // convert_map_to_range_list function below.
      for (auto it : vector_partitioner->import_indices())
        for (unsigned int i = it.first; i < it.second; ++i)
          touched_first_by[i / chunk_size_zero_vector] = n_partitions;
      std::map<unsigned int, std::vector<unsigned int>> chunk_must_do_pre;
      for (unsigned int i = 0; i < touched_first_by.size(); ++i)
        chunk_must_do_pre[touched_first_by[i]].push_back(i);
      convert_map_to_range_list(n_partitions + 1,
                                chunk_must_do_pre,
                                cell_loop_pre_list_index,
                                cell_loop_pre_list,
                                vector_partitioner->locally_owned_size());

      touched_last_by.resize((vector_partitioner->locally_owned_size() +
                              chunk_size_zero_vector - 1) /
                             chunk_size_zero_vector);

      // set the indices which were not touched by the cell loop (i.e.,
      // constrained indices) to the last valid partition index. Since
      // partition_row_index contains one extra slot for ghosted faces (which
      // are not part of the cell/face loops), we use the second to last entry
      // in the partition list.
      for (auto &index : touched_last_by)
        if (index == numbers::invalid_unsigned_int)
          index =
            task_info
              .partition_row_index[task_info.partition_row_index.size() - 2] -
            1;
      for (auto it : vector_partitioner->import_indices())
        for (unsigned int i = it.first; i < it.second; ++i)
          touched_last_by[i / chunk_size_zero_vector] = n_partitions;
      std::map<unsigned int, std::vector<unsigned int>> chunk_must_do_post;
      for (unsigned int i = 0; i < touched_last_by.size(); ++i)
        chunk_must_do_post[touched_last_by[i]].push_back(i);
      convert_map_to_range_list(n_partitions + 1,
                                chunk_must_do_post,
                                cell_loop_post_list_index,
                                cell_loop_post_list,
                                vector_partitioner->locally_owned_size());
    }



    namespace internal
    {
      // rudimentary version of a vector that keeps entries always ordered
      class ordered_vector : public std::vector<types::global_dof_index>
      {
      public:
        ordered_vector()
        {
          reserve(2000);
        }

        void
        reserve(const std::size_t size)
        {
          if (size > 0)
            this->std::vector<types::global_dof_index>::reserve(size);
        }


        // insert a given entry. dat is a pointer within this vector (the user
        // needs to make sure that it really stays there)
        void
        insert(const unsigned int                              entry,
               std::vector<types::global_dof_index>::iterator &dat)
        {
          AssertIndexRange(static_cast<std::size_t>(dat - begin()), size() + 1);
          AssertIndexRange(static_cast<std::size_t>(end() - dat), size() + 1);
          AssertIndexRange(size(), capacity());
          while (dat != end() && *dat < entry)
            ++dat;

          if (dat == end())
            {
              push_back(entry);
              dat = end();
            }
          else if (*dat > entry)
            {
              dat =
                this->std::vector<types::global_dof_index>::insert(dat, entry);
              ++dat;
            }
          else
            ++dat;
        }
      };
    } // namespace internal



    template <typename StreamType>
    void
    DoFInfo::print_memory_consumption(StreamType     &out,
                                      const TaskInfo &task_info) const
    {
      out << "       Memory row starts indices:    ";
      task_info.print_memory_statistics(out,
                                        (row_starts.capacity() *
                                         sizeof(*row_starts.begin())));
      out << "       Memory dof indices:           ";
      task_info.print_memory_statistics(
        out, MemoryConsumption::memory_consumption(dof_indices));
      out << "       Memory constraint indicators: ";
      task_info.print_memory_statistics(
        out, MemoryConsumption::memory_consumption(constraint_indicator));
      out << "       Memory plain indices:         ";
      task_info.print_memory_statistics(
        out,
        MemoryConsumption::memory_consumption(row_starts_plain_indices) +
          MemoryConsumption::memory_consumption(plain_dof_indices));
      out << "       Memory vector partitioner:    ";
      task_info.print_memory_statistics(
        out, MemoryConsumption::memory_consumption(*vector_partitioner));
    }



    template <typename Number>
    void
    DoFInfo::print(const std::vector<Number>       &constraint_pool_data,
                   const std::vector<unsigned int> &constraint_pool_row_index,
                   std::ostream                    &out) const
    {
      const unsigned int n_rows = row_starts.size() - 1;
      for (unsigned int row = 0; row < n_rows; ++row)
        {
          if (row_starts[row].first == row_starts[row + 1].first)
            continue;
          out << "Entries row " << row << ": ";
          const unsigned int *glob_indices =
                               &dof_indices[row_starts[row].first],
                             *end_row = &dof_indices[row_starts[row + 1].first];
          unsigned int index          = 0;
          const std::pair<unsigned short, unsigned short>
            *con_it  = &constraint_indicator[row_starts[row].second],
            *end_con = &constraint_indicator[row_starts[row + 1].second];
          for (; con_it != end_con; ++con_it)
            {
              for (unsigned int j = 0; j < con_it->first; ++j, ++index)
                {
                  Assert(glob_indices + index != end_row, ExcInternalError());
                  out << glob_indices[index] << ' ';
                }

              out << "[ ";
              for (unsigned int k = constraint_pool_row_index[con_it->second];
                   k < constraint_pool_row_index[con_it->second + 1];
                   k++, index++)
                {
                  Assert(glob_indices + index != end_row, ExcInternalError());
                  out << glob_indices[index] << '/' << constraint_pool_data[k]
                      << ' ';
                }
              out << "] ";
            }
          glob_indices += index;
          for (; glob_indices != end_row; ++glob_indices)
            out << *glob_indices << ' ';
          out << std::endl;
        }
    }


  } // namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
