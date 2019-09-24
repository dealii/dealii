// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2019 by the deal.II authors
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

#ifndef dealii_matrix_free_dof_info_templates_h
#define dealii_matrix_free_dof_info_templates_h


#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/mapping_info.h>
#include <deal.II/matrix_free/task_info.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    struct ConstraintComparator
    {
      bool
      operator()(const std::pair<types::global_dof_index, double> &p1,
                 const std::pair<types::global_dof_index, double> &p2) const
      {
        return p1.second < p2.second;
      }
    };

    /**
     * A struct that takes entries describing a constraint and puts them into
     * a sorted list where duplicates are filtered out
     */
    template <typename Number>
    struct ConstraintValues
    {
      ConstraintValues();

      /**
       * This function inserts some constrained entries to the collection of
       * all values. It stores the (reordered) numbering of the dofs
       * (according to the ordering that matches with the function) in
       * new_indices, and returns the storage position the double array for
       * access later on.
       */
      template <typename number2>
      unsigned short
      insert_entries(
        const std::vector<std::pair<types::global_dof_index, number2>>
          &entries);

      std::vector<std::pair<types::global_dof_index, double>>
                                           constraint_entries;
      std::vector<types::global_dof_index> constraint_indices;

      std::pair<std::vector<Number>, types::global_dof_index> next_constraint;
      std::map<std::vector<Number>,
               types::global_dof_index,
               FPArrayComparator<double>>
        constraints;
    };


    template <typename Number>
    ConstraintValues<Number>::ConstraintValues()
      : constraints(FPArrayComparator<double>(1.))
    {}

    template <typename Number>
    template <typename number2>
    unsigned short
    ConstraintValues<Number>::insert_entries(
      const std::vector<std::pair<types::global_dof_index, number2>> &entries)
    {
      next_constraint.first.resize(entries.size());
      if (entries.size() > 0)
        {
          constraint_indices.resize(entries.size());
          // Use assign so that values for nonmatching Number / number2 are
          // converted:
          constraint_entries.assign(entries.begin(), entries.end());
          std::sort(constraint_entries.begin(),
                    constraint_entries.end(),
                    ConstraintComparator());
          for (types::global_dof_index j = 0; j < constraint_entries.size();
               j++)
            {
              // copy the indices of the constraint entries after sorting.
              constraint_indices[j] = constraint_entries[j].first;

              // one_constraint takes the weights of the constraint
              next_constraint.first[j] = constraint_entries[j].second;
            }
        }
      next_constraint.second = constraints.size();

      // check whether or not constraint is already in pool. the initial
      // implementation computed a hash value based on the truncated array (to
      // given accuracy around 1e-13) in order to easily detect different
      // arrays and then made a fine-grained check when the hash values were
      // equal. this was quite lengthy and now we use a std::map with a
      // user-defined comparator to compare floating point arrays to a
      // tolerance 1e-13.
      std::pair<typename std::map<std::vector<double>,
                                  types::global_dof_index,
                                  FPArrayComparator<double>>::iterator,
                bool>
        it = constraints.insert(next_constraint);

      types::global_dof_index insert_position = numbers::invalid_dof_index;
      if (it.second == false)
        insert_position = it.first->second;
      else
        insert_position = next_constraint.second;

      // we want to store the result as a short variable, so we have to make
      // sure that the result does not exceed the limits when casting.
      Assert(insert_position < (1 << (8 * sizeof(unsigned short))),
             ExcInternalError());
      return static_cast<unsigned short>(insert_position);
    }



    // ----------------- actual DoFInfo functions -----------------------------

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
      dimension                  = 2;
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
            row_starts[ib].second != row_starts[ib + n_fe_components].second;

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


    template <typename number>
    void
    DoFInfo ::read_dof_indices(
      const std::vector<types::global_dof_index> &local_indices,
      const std::vector<unsigned int> &           lexicographic_inv,
      const AffineConstraints<number> &           constraints,
      const unsigned int                          cell_number,
      ConstraintValues<double> &                  constraint_values,
      bool &                                      cell_at_subdomain_boundary)
    {
      Assert(vector_partitioner.get() != nullptr, ExcInternalError());
      const unsigned int n_mpi_procs = vector_partitioner->n_mpi_processes();
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
      for (unsigned int comp = 0; comp < n_components; ++comp)
        {
          std::pair<unsigned short, unsigned short> constraint_iterator(0, 0);
          for (unsigned int i = component_dof_indices_offset[fe_index][comp];
               i < component_dof_indices_offset[fe_index][comp + 1];
               i++)
            {
              types::global_dof_index current_dof =
                local_indices[lexicographic_inv[i]];
              const auto *entries_ptr =
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
                  const auto &                  entries   = *entries_ptr;
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
                          if (n_mpi_procs > 1 &&
                              (constraint_indices[j] < first_owned ||
                               constraint_indices[j] >= last_owned))
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
                  if (n_mpi_procs > 1 &&
                      (current_dof < first_owned || current_dof >= last_owned))
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
            (row_starts[(cell_number + 1) * n_components].second >
             row_starts[cell_number * n_components].second);
          if (cell_has_constraints == true)
            {
              for (unsigned int i = 0; i < dofs_this_cell; ++i)
                {
                  types::global_dof_index current_dof =
                    local_indices[lexicographic_inv[i]];
                  if (n_mpi_procs > 1 &&
                      (current_dof < first_owned || current_dof >= last_owned))
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



    void
    DoFInfo ::assign_ghosts(const std::vector<unsigned int> &boundary_cells)
    {
      Assert(boundary_cells.size() < row_starts.size(), ExcInternalError());

      // sort ghost dofs and compress out duplicates
      const unsigned int n_owned  = (vector_partitioner->local_range().second -
                                    vector_partitioner->local_range().first);
      const std::size_t  n_ghosts = ghost_dofs.size();
#ifdef DEBUG
      for (const auto dof_index : dof_indices)
        AssertIndexRange(dof_index, n_owned + n_ghosts);
#endif

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
          for (std::size_t i = 1; i < n_ghosts; i++)
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
                  if (row_starts[boundary_cells[i] * n_components].second !=
                      row_starts[(boundary_cells[i] + 1) * n_components].second)
                    {
                      unsigned int *data_ptr =
                        plain_dof_indices.data() +
                        row_starts_plain_indices[boundary_cells[i]];
                      const unsigned int fe_index =
                        (cell_active_fe_index.size() == 0 ||
                         dofs_per_cell.size() == 1) ?
                          0 :
                          cell_active_fe_index[i];
                      AssertIndexRange(fe_index, dofs_per_cell.size());
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
    }



    void
    DoFInfo ::reorder_cells(
      const TaskInfo &                  task_info,
      const std::vector<unsigned int> & renumbering,
      const std::vector<unsigned int> & constraint_pool_row_index,
      const std::vector<unsigned char> &irregular_cells)
    {
      (void)constraint_pool_row_index;

      // first reorder the active fe index.
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
                  row_starts[cell_no].second !=
                    row_starts[cell_no + n_components].second)
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

#ifdef DEBUG
      // sanity check 1: all indices should be smaller than the number of dofs
      // locally owned plus the number of ghosts
      const unsigned int index_range =
        (vector_partitioner->local_range().second -
         vector_partitioner->local_range().first) +
        vector_partitioner->ghost_indices().n_elements();
      for (const auto dof_index : dof_indices)
        AssertIndexRange(dof_index, index_range);

      // sanity check 2: for the constraint indicators, the first index should
      // be smaller than the number of indices in the row, and the second
      // index should be smaller than the number of constraints in the
      // constraint pool.
      for (unsigned int row = 0; row < task_info.cell_partition_data.back();
           ++row)
        {
          const unsigned int row_length_ind =
            row_starts[(row * vectorization_length + 1) * n_components].first -
            row_starts[row * vectorization_length * n_components].first;
          AssertIndexRange(
            row_starts[(row * vectorization_length + 1) * n_components].second,
            constraint_indicator.size() + 1);
          const std::pair<unsigned short, unsigned short> *
            con_it =
             constraint_indicator.data() +
             row_starts[row * vectorization_length * n_components].second,
           *end_con =
             constraint_indicator.data() +
             row_starts[(row * vectorization_length + 1) * n_components].second;
          for (; con_it != end_con; ++con_it)
            {
              AssertIndexRange(con_it->first, row_length_ind + 1);
              AssertIndexRange(con_it->second,
                               constraint_pool_row_index.size() - 1);
            }
        }

      // sanity check 3: check the number of cells once again
      unsigned int n_active_cells = 0;
      for (unsigned int c = 0; c < *(task_info.cell_partition_data.end() - 2);
           ++c)
        if (irregular_cells[c] > 0)
          n_active_cells += irregular_cells[c];
        else
          n_active_cells += vectorization_length;
      AssertDimension(n_active_cells, task_info.n_active_cells);
#endif

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
              bool indices_are_contiguous = true;
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
                for (unsigned int k = 0; k < ndofs; ++k)
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
                    dof_indices_contiguous
                      [dof_access_cell][i * vectorization_length + j] =
                        this->dof_indices[row_starts[(i * vectorization_length +
                                                      j) *
                                                     n_components]
                                            .first];
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
              else
                {
                  int                 indices_are_interleaved_and_mixed = 2;
                  const unsigned int *dof_indices =
                    &this->dof_indices[row_starts[i * vectorization_length *
                                                  n_components]
                                         .first];
                  for (unsigned int j = 0; j < n_comp; ++j)
                    offsets[j] =
                      dof_indices[j * ndofs + 1] - dof_indices[j * ndofs];
                  for (unsigned int k = 0; k < ndofs; ++k)
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

                      // do not use interleaved storage if two vectorized
                      // components point to the same field (scatter not
                      // possible)
                      for (unsigned int k = 0; k < ndofs; ++k)
                        for (unsigned int l = 0; l < n_comp; ++l)
                          for (unsigned int j = l + 1; j < n_comp; ++j)
                            if (dof_indices[j * ndofs + k] ==
                                dof_indices[l * ndofs + k])
                              {
                                index_storage_variants[dof_access_cell][i] =
                                  IndexStorageVariants::full;
                                break;
                              }
                    }
                }
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
              this->dof_indices_interleaved.size());
            AssertIndexRange(
              row_starts[i * vectorization_length * n_components].first +
                ndofs * vectorization_length,
              this->dof_indices_interleaved.size() + 1);
            for (unsigned int k = 0; k < ndofs; ++k)
              for (unsigned int j = 0; j < vectorization_length; ++j)
                interleaved_dof_indices[k * vectorization_length + j] =
                  dof_indices[j * ndofs + k];
          }
    }



    template <int length>
    void
    DoFInfo::compute_face_index_compression(
      const std::vector<FaceToCellTopology<length>> &faces)
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

      // all interior faces come before the boundary faces
      unsigned int n_exterior_faces = 0;
      for (; n_exterior_faces < faces.size(); ++n_exterior_faces)
        if (faces[n_exterior_faces].cells_exterior[0] ==
            numbers::invalid_unsigned_int)
          break;
      index_storage_variants[dof_access_face_exterior].resize(
        n_exterior_faces, IndexStorageVariants::full);
      dof_indices_contiguous[dof_access_face_exterior].resize(
        n_exterior_faces * length, numbers::invalid_unsigned_int);
      dof_indices_interleave_strides[dof_access_face_exterior].resize(
        faces.size() * length, numbers::invalid_unsigned_int);
      n_vectorization_lanes_filled[dof_access_face_exterior].resize(
        n_exterior_faces);

      for (unsigned int face = 0; face < faces.size(); ++face)
        {
          auto face_computation = [&](const DoFAccessIndex face_index,
                                      const unsigned int * cell_indices_face) {
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
          if (face < n_exterior_faces)
            face_computation(dof_access_face_exterior,
                             faces[face].cells_exterior);
        }
    }



    template <int length>
    void
    DoFInfo::compute_vector_zero_access_pattern(
      const TaskInfo &                               task_info,
      const std::vector<FaceToCellTopology<length>> &faces)
    {
      // compute a list that tells us the first time a degree of freedom is
      // touched by a cell
      AssertDimension(length, vectorization_length);
      const unsigned int n_components = start_components.back();
      const unsigned int n_dofs       = vector_partitioner->local_size() +
                                  vector_partitioner->n_ghost_indices();
      std::vector<unsigned int> touched_first_by(
        (n_dofs + chunk_size_zero_vector - 1) / chunk_size_zero_vector,
        numbers::invalid_unsigned_int);
      std::vector<unsigned int> touched_last_by(
        (n_dofs + chunk_size_zero_vector - 1) / chunk_size_zero_vector,
        numbers::invalid_unsigned_int);
      for (unsigned int part = 0;
           part < task_info.partition_row_index.size() - 2;
           ++part)
        for (unsigned int chunk = task_info.partition_row_index[part];
             chunk < task_info.partition_row_index[part + 1];
             ++chunk)
          {
            for (unsigned int cell = task_info.cell_partition_data[chunk];
                 cell < task_info.cell_partition_data[chunk + 1];
                 ++cell)
              {
                for (unsigned int it =
                       row_starts[cell * vectorization_length * n_components]
                         .first;
                     it != row_starts[(cell + 1) * vectorization_length *
                                      n_components]
                             .first;
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
            if (faces.size() > 0)
              for (unsigned int face = task_info.face_partition_data[chunk];
                   face < task_info.face_partition_data[chunk + 1];
                   ++face)
                for (unsigned int v = 0;
                     v < length && faces[face].cells_exterior[v] !=
                                     numbers::invalid_unsigned_int;
                     ++v)
                  {
                    const unsigned int cell = faces[face].cells_exterior[v];
                    for (unsigned int it =
                           row_starts[cell * n_components].first;
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
                                vector_partitioner->local_size());

      // the other two operations only work on the local range (without
      // ghosts), so we skip the latter parts of the vector now
      touched_first_by.resize(
        (vector_partitioner->local_size() + chunk_size_zero_vector - 1) /
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
                                vector_partitioner->local_size());

      touched_last_by.resize(
        (vector_partitioner->local_size() + chunk_size_zero_vector - 1) /
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
                                vector_partitioner->local_size());
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



      // We construct the connectivity graph in parallel. we use one lock for
      // 256 degrees of freedom to keep the number of locks down to a
      // reasonable level and reduce the cost of locking to some extent.
      static constexpr unsigned int bucket_size_threading = 256;

      void
      compute_row_lengths(const unsigned int           begin,
                          const unsigned int           end,
                          const DoFInfo &              dof_info,
                          std::vector<Threads::Mutex> &mutexes,
                          std::vector<unsigned int> &  row_lengths)
      {
        std::vector<unsigned int> scratch;
        const unsigned int n_components = dof_info.start_components.back();
        for (unsigned int block = begin; block < end; ++block)
          {
            scratch.clear();
            scratch.insert(
              scratch.end(),
              dof_info.dof_indices.data() +
                dof_info.row_starts[block * n_components].first,
              dof_info.dof_indices.data() +
                dof_info.row_starts[(block + 1) * n_components].first);
            std::sort(scratch.begin(), scratch.end());
            std::vector<unsigned int>::const_iterator end_unique =
              std::unique(scratch.begin(), scratch.end());
            std::vector<unsigned int>::const_iterator it = scratch.begin();
            while (it != end_unique)
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
                    row_lengths[*it]++;
                  }
              }
          }
      }

      void
      fill_connectivity_dofs(const unsigned int               begin,
                             const unsigned int               end,
                             const DoFInfo &                  dof_info,
                             const std::vector<unsigned int> &row_lengths,
                             std::vector<Threads::Mutex> &    mutexes,
                             dealii::SparsityPattern &        connectivity_dof)
      {
        std::vector<unsigned int> scratch;
        const unsigned int n_components = dof_info.start_components.back();
        for (unsigned int block = begin; block < end; ++block)
          {
            scratch.clear();
            scratch.insert(
              scratch.end(),
              dof_info.dof_indices.data() +
                dof_info.row_starts[block * n_components].first,
              dof_info.dof_indices.data() +
                dof_info.row_starts[(block + 1) * n_components].first);
            std::sort(scratch.begin(), scratch.end());
            std::vector<unsigned int>::const_iterator end_unique =
              std::unique(scratch.begin(), scratch.end());
            std::vector<unsigned int>::const_iterator it = scratch.begin();
            while (it != end_unique)
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
                        const DoFInfo &                  dof_info,
                        const std::vector<unsigned int> &renumbering,
                        const dealii::SparsityPattern &  connectivity_dof,
                        DynamicSparsityPattern &         connectivity)
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
      const TaskInfo &                 task_info,
      const std::vector<unsigned int> &renumbering,
      DynamicSparsityPattern &         connectivity) const
    {
      unsigned int n_rows = (vector_partitioner->local_range().second -
                             vector_partitioner->local_range().first) +
                            vector_partitioner->ghost_indices().n_elements();

      // Avoid square sparsity patterns that allocate the diagonal entry
      if (n_rows == task_info.n_active_cells)
        ++n_rows;

      // first determine row lengths
      std::vector<unsigned int>   row_lengths(n_rows);
      std::vector<Threads::Mutex> mutexes(
        n_rows / internal::bucket_size_threading + 1);
      parallel::apply_to_subranges(0,
                                   task_info.n_active_cells,
                                   [this,
                                    &mutexes,
                                    &row_lengths](const unsigned int begin,
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
      parallel::apply_to_subranges(
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
      parallel::apply_to_subranges(
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
    DoFInfo ::compute_dof_renumbering(
      std::vector<types::global_dof_index> &renumbering)
    {
      const unsigned int local_size = vector_partitioner->local_size();
      renumbering.resize(0);
      renumbering.resize(local_size, numbers::invalid_dof_index);

      types::global_dof_index counter      = 0;
      const unsigned int      n_components = start_components.back();
      const unsigned int      n_macro_cells =
        n_vectorization_lanes_filled[dof_access_cell].size();
      Assert(n_macro_cells <=
               (row_starts.size() - 1) / vectorization_length / n_components,
             ExcInternalError());
      for (unsigned int cell_no = 0; cell_no < n_macro_cells; ++cell_no)
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
                  if (dof_ind[j * ndofs + i] < local_size)
                    if (renumbering[dof_ind[j * ndofs + i]] ==
                        numbers::invalid_dof_index)
                      renumbering[dof_ind[j * ndofs + i]] = counter++;
            }
        }

      AssertIndexRange(counter, local_size + 1);
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
      memory +=
        (row_starts.capacity() * sizeof(std::pair<unsigned int, unsigned int>));
      memory += MemoryConsumption::memory_consumption(dof_indices);
      memory += MemoryConsumption::memory_consumption(row_starts_plain_indices);
      memory += MemoryConsumption::memory_consumption(plain_dof_indices);
      memory += MemoryConsumption::memory_consumption(constraint_indicator);
      memory += MemoryConsumption::memory_consumption(*vector_partitioner);
      return memory;
    }



    template <typename StreamType>
    void
    DoFInfo::print_memory_consumption(StreamType &    out,
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
    DoFInfo::print(const std::vector<Number> &      constraint_pool_data,
                   const std::vector<unsigned int> &constraint_pool_row_index,
                   std::ostream &                   out) const
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
                  out << glob_indices[index] << " ";
                }

              out << "[ ";
              for (unsigned int k = constraint_pool_row_index[con_it->second];
                   k < constraint_pool_row_index[con_it->second + 1];
                   k++, index++)
                {
                  Assert(glob_indices + index != end_row, ExcInternalError());
                  out << glob_indices[index] << "/" << constraint_pool_data[k]
                      << " ";
                }
              out << "] ";
            }
          glob_indices += index;
          for (; glob_indices != end_row; ++glob_indices)
            out << *glob_indices << " ";
          out << std::endl;
        }
    }


  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
