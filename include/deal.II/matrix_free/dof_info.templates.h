// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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

#ifndef dealii_matrix_free_dof_info_templates_h
#define dealii_matrix_free_dof_info_templates_h


#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/helper_functions.h>
#include <deal.II/matrix_free/mapping_info.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {

    struct ConstraintComparator
    {
      bool operator()(const std::pair<types::global_dof_index,double> &p1,
                      const std::pair<types::global_dof_index,double> &p2) const
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
      unsigned short
      insert_entries (const std::vector<std::pair<types::global_dof_index,double> > &entries);

      std::vector<std::pair<types::global_dof_index, double> > constraint_entries;
      std::vector<types::global_dof_index> constraint_indices;

      std::pair<std::vector<Number>, types::global_dof_index> next_constraint;
      std::map<std::vector<Number>, types::global_dof_index, FPArrayComparator<double> > constraints;
    };


    template <typename Number>
    ConstraintValues<Number>::ConstraintValues ()
      :
      constraints(FPArrayComparator<double>(1.))
    {}

    template <typename Number>
    unsigned short
    ConstraintValues<Number>::
    insert_entries (const std::vector<std::pair<types::global_dof_index,double> > &entries)
    {
      next_constraint.first.resize(entries.size());
      if (entries.size() > 0)
        {
          constraint_indices.resize(entries.size());
          constraint_entries = entries;
          std::sort(constraint_entries.begin(), constraint_entries.end(),
                    ConstraintComparator());
          for (types::global_dof_index j=0; j<constraint_entries.size(); j++)
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
      std::pair<typename std::map<std::vector<double>, types::global_dof_index,
          FPArrayComparator<double> >::iterator,
          bool> it = constraints.insert(next_constraint);

      types::global_dof_index insert_position = numbers::invalid_dof_index;
      if (it.second == false)
        insert_position = it.first->second;
      else
        insert_position = next_constraint.second;

      // we want to store the result as a short variable, so we have to make
      // sure that the result does not exceed the limits when casting.
      Assert(insert_position < (1<<(8*sizeof(unsigned short))),
             ExcInternalError());
      return static_cast<unsigned short>(insert_position);
    }



    // ----------------- actual DoFInfo functions -----------------------------

    DoFInfo::DoFInfo ()
    {
      clear();
    }



    void
    DoFInfo::clear ()
    {
      row_starts.clear();
      dof_indices.clear();
      constraint_indicator.clear();
      vector_partitioner.reset();
      ghost_dofs.clear();
      dofs_per_cell.clear();
      dofs_per_face.clear();
      dimension = 2;
      n_components = 0;
      row_starts_plain_indices.clear();
      plain_dof_indices.clear();
      store_plain_indices = false;
      cell_active_fe_index.clear();
      max_fe_index = 0;
      fe_index_conversion.clear();
    }



    void
    DoFInfo::read_dof_indices (const std::vector<types::global_dof_index> &local_indices,
                               const std::vector<unsigned int> &lexicographic_inv,
                               const ConstraintMatrix          &constraints,
                               const unsigned int               cell_number,
                               ConstraintValues<double> &constraint_values,
                               bool                            &cell_at_boundary)
    {
      Assert (vector_partitioner.get() != nullptr, ExcInternalError());
      const unsigned int n_mpi_procs = vector_partitioner->n_mpi_processes();
      const types::global_dof_index first_owned = vector_partitioner->local_range().first;
      const types::global_dof_index last_owned  = vector_partitioner->local_range().second;
      Assert (last_owned-first_owned < std::numeric_limits<unsigned int>::max(),
              ExcMessage("The size local range of owned indices must not "
                         "exceed the size of unsigned int"));
      const unsigned int n_owned     = last_owned - first_owned;
      std::pair<unsigned short,unsigned short> constraint_iterator (0,0);

      unsigned int dofs_this_cell = (cell_active_fe_index.empty()) ?
                                    dofs_per_cell[0] : dofs_per_cell[cell_active_fe_index[cell_number]];
      for (unsigned int i=0; i<dofs_this_cell; i++)
        {
          types::global_dof_index current_dof =
            local_indices[lexicographic_inv[i]];
          const std::vector<std::pair<types::global_dof_index,double> >
          *entries_ptr =
            constraints.get_constraint_entries(current_dof);

          // dof is constrained
          if (entries_ptr != nullptr)
            {
              // in case we want to access plain indices, we need to know
              // about the location of constrained indices as well (all the
              // other indices are collected by the cases below)
              if (current_dof < first_owned || current_dof >= last_owned)
                {
                  ghost_dofs.push_back (current_dof);
                  cell_at_boundary = true;
                }

              // check whether this dof is identity constrained to another
              // dof. then we can simply insert that dof and there is no need
              // to actually resolve the constraint entries
              const std::vector<std::pair<types::global_dof_index,double> >
              &entries = *entries_ptr;
              const types::global_dof_index n_entries = entries.size();
              if (n_entries == 1 && std::fabs(entries[0].second-1.)<1e-14)
                {
                  current_dof = entries[0].first;
                  goto no_constraint;
                }

              // append a new index to the indicators
              constraint_indicator.push_back (constraint_iterator);
              constraint_indicator.back().second =
                constraint_values.insert_entries (entries);

              // reset constraint iterator for next round
              constraint_iterator.first = 0;

              // add the local_to_global indices computed in the
              // insert_entries function. transform the index to local index
              // space or mark it as ghost if necessary
              if (n_entries > 0)
                {
                  const std::vector<types::global_dof_index> &constraint_indices =
                    constraint_values.constraint_indices;
                  for (unsigned int j=0; j<n_entries; ++j)
                    {
                      if (n_mpi_procs > 1 &&
                          (constraint_indices[j] < first_owned ||
                           constraint_indices[j] >= last_owned))
                        {
                          dof_indices.push_back (n_owned + ghost_dofs.size());

                          // collect ghosts so that we can later construct an
                          // IndexSet for them. also store whether the current
                          // cell is on the boundary
                          ghost_dofs.push_back(constraint_indices[j]);
                          cell_at_boundary = true;
                        }
                      else
                        // not ghost, so transform to the local index space
                        // directly
                        dof_indices.push_back
                        (static_cast<unsigned int>(constraint_indices[j] -
                                                   first_owned));
                    }
                }
            }
          else
            {
no_constraint:
              // Not constrained, we simply have to add the local index to the
              // indices_local_to_global list and increment constraint
              // iterator. transform to local index space/mark as ghost
              if (n_mpi_procs > 1 &&
                  (current_dof < first_owned ||
                   current_dof >= last_owned))
                {
                  ghost_dofs.push_back(current_dof);
                  current_dof = n_owned + ghost_dofs.size()-1;
                  cell_at_boundary = true;
                }
              else
                current_dof -= first_owned;

              dof_indices.push_back (static_cast<unsigned int>(current_dof));

              // make sure constraint_iterator.first is always within the
              // bounds of unsigned short
              Assert (constraint_iterator.first <
                      (1<<(8*sizeof(unsigned short)))-1,
                      ExcInternalError());
              constraint_iterator.first++;
            }
        }
      row_starts[cell_number+1][0] = dof_indices.size();
      row_starts[cell_number+1][1] = constraint_indicator.size();
      row_starts[cell_number+1][2] = 0;

      // now to the plain indices: in case we have constraints on this cell,
      // store the indices without the constraints resolve once again
      if (store_plain_indices == true)
        {
          if (cell_number == 0)
            row_starts_plain_indices.resize (row_starts.size());
          row_starts_plain_indices[cell_number] = plain_dof_indices.size();
          bool cell_has_constraints = (row_starts[cell_number+1][1] >
                                       row_starts[cell_number][1]);
          if (cell_has_constraints == true)
            {
              for (unsigned int i=0; i<dofs_this_cell; ++i)
                {
                  types::global_dof_index current_dof =
                    local_indices[lexicographic_inv[i]];
                  if (n_mpi_procs > 1 &&
                      (current_dof < first_owned ||
                       current_dof >= last_owned))
                    {
                      ghost_dofs.push_back(current_dof);
                      current_dof = n_owned + ghost_dofs.size()-1;
                      cell_at_boundary = true;
                    }
                  else
                    current_dof -= first_owned;
                  plain_dof_indices.push_back (static_cast<unsigned int>
                                               (current_dof));
                }
            }
        }
    }



    void
    DoFInfo::assign_ghosts (const std::vector<unsigned int> &boundary_cells)
    {
      Assert (boundary_cells.size() < row_starts.size(), ExcInternalError());

      // sort ghost dofs and compress out duplicates
      const unsigned int n_owned  = (vector_partitioner->local_range().second-
                                     vector_partitioner->local_range().first);
      const std::size_t n_ghosts = ghost_dofs.size();
#ifdef DEBUG
      for (std::vector<unsigned int>::iterator dof = dof_indices.begin();
           dof!=dof_indices.end(); ++dof)
        AssertIndexRange (*dof, n_owned+n_ghosts);
#endif

      std::vector<unsigned int> ghost_numbering (n_ghosts);
      IndexSet ghost_indices (vector_partitioner->size());
      if (n_ghosts > 0)
        {
          unsigned int n_unique_ghosts = 0;
          // since we need to go back to the local_to_global indices and
          // replace the temporary numbering of ghosts by the real number in
          // the index set, we need to store these values
          std::vector<std::pair<types::global_dof_index,unsigned int> > ghost_origin(n_ghosts);
          for (std::size_t i=0; i<n_ghosts; ++i)
            {
              ghost_origin[i].first = ghost_dofs[i];
              ghost_origin[i].second = i;
            }
          std::sort (ghost_origin.begin(), ghost_origin.end());

          types::global_dof_index last_contiguous_start = ghost_origin[0].first;
          ghost_numbering[ghost_origin[0].second] = 0;
          for (std::size_t i=1; i<n_ghosts; i++)
            {
              if (ghost_origin[i].first > ghost_origin[i-1].first+1)
                {
                  ghost_indices.add_range (last_contiguous_start,
                                           ghost_origin[i-1].first+1);
                  last_contiguous_start = ghost_origin[i].first;
                }
              if (ghost_origin[i].first>ghost_origin[i-1].first)
                ++n_unique_ghosts;
              ghost_numbering[ghost_origin[i].second] = n_unique_ghosts;
            }
          ++n_unique_ghosts;
          ghost_indices.add_range (last_contiguous_start,
                                   ghost_origin.back().first+1);
          ghost_indices.compress();

          // make sure that we got the correct local numbering of the ghost
          // dofs. the ghost index set should store the same number
          {
            AssertDimension (n_unique_ghosts, ghost_indices.n_elements());
            for (std::size_t i=0; i<n_ghosts; ++i)
              Assert (ghost_numbering[i] ==
                      ghost_indices.index_within_set(ghost_dofs[i]),
                      ExcInternalError());
          }

          // apply correct numbering for ghost indices: We previously just
          // enumerated them according to their appearance in the
          // local_to_global structure. Above, we derived a relation between
          // this enumeration and the actual number
          const unsigned int n_boundary_cells = boundary_cells.size();
          for (unsigned int i=0; i<n_boundary_cells; ++i)
            {
              unsigned int *data_ptr = const_cast<unsigned int *> (begin_indices(boundary_cells[i]));

              const unsigned int *row_end = end_indices(boundary_cells[i]);
              for ( ; data_ptr != row_end; ++data_ptr)
                *data_ptr = ((*data_ptr < n_owned)
                             ?
                             *data_ptr
                             :
                             n_owned +
                             ghost_numbering[*data_ptr - n_owned]);

              // now the same procedure for plain indices
              if (store_plain_indices == true)
                {
                  if (row_length_indicators(boundary_cells[i]) > 0)
                    {
                      unsigned int *data_ptr = const_cast<unsigned int *> (begin_indices_plain(boundary_cells[i]));
                      const unsigned int *row_end = end_indices_plain(boundary_cells[i]);
                      for ( ; data_ptr != row_end; ++data_ptr)
                        *data_ptr = ((*data_ptr < n_owned)
                                     ?
                                     *data_ptr
                                     :
                                     n_owned +
                                     ghost_numbering[*data_ptr - n_owned]);
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
      vec_part->set_ghost_indices (ghost_indices);
    }



    void
    DoFInfo::reorder_cells (const TaskInfo                   &task_info,
                            const std::vector<unsigned int>  &renumbering,
                            const std::vector<unsigned int>  &constraint_pool_row_index,
                            const std::vector<unsigned char> &irregular_cells,
                            const unsigned int                vectorization_length)
    {
      // first reorder the active fe index.
      if (cell_active_fe_index.size() > 0)
        {
          std::vector<unsigned int> new_active_fe_index;
          new_active_fe_index.reserve (task_info.cell_partition_data.back());
          std::vector<unsigned int> fe_indices(vectorization_length);
          unsigned int position_cell = 0;
          for (unsigned int cell=0; cell<task_info.cell_partition_data.back(); ++cell)
            {
              const unsigned int n_comp = (irregular_cells[cell] > 0 ?
                                           irregular_cells[cell] : vectorization_length);
              for (unsigned int j=0; j<n_comp; ++j)
                fe_indices[j]=cell_active_fe_index[renumbering[position_cell+j]];

              // by construction, all cells should have the same fe index.
              for (unsigned int j=1; j<n_comp; ++j)
                Assert (fe_indices[j] == fe_indices[0], ExcInternalError());

              new_active_fe_index.push_back(fe_indices[0]);
              position_cell += n_comp;
            }
          std::swap (new_active_fe_index, cell_active_fe_index);
        }

      std::vector<std::array<unsigned int, 3> > new_row_starts;
      std::vector<unsigned int> new_dof_indices;
      std::vector<std::pair<unsigned short,unsigned short> >
      new_constraint_indicator;
      std::vector<unsigned int> new_plain_indices, new_rowstart_plain;
      unsigned int position_cell = 0;
      new_row_starts.resize(task_info.cell_partition_data.back()+1);
      new_dof_indices.reserve (dof_indices.size());
      new_constraint_indicator.reserve (constraint_indicator.size());
      if (store_plain_indices == true)
        {
          new_rowstart_plain.resize (task_info.cell_partition_data.back()+1,
                                     numbers::invalid_unsigned_int);
          new_plain_indices.reserve (plain_dof_indices.size());
        }

      // copy the indices and the constraint indicators to the new data field:
      // Store the indices in a way so that adjacent data fields in local
      // vectors are adjacent, i.e., first dof index 0 for all vectors, then
      // dof index 1 for all vectors, and so on. This involves some extra
      // resorting.
      std::vector<const unsigned int *> glob_indices (vectorization_length);
      std::vector<const unsigned int *> plain_glob_indices (vectorization_length);
      std::vector<const std::pair<unsigned short,unsigned short>*>
      constr_ind(vectorization_length), constr_end(vectorization_length);
      std::vector<unsigned int> index(vectorization_length);
      for (unsigned int i=0; i<task_info.cell_partition_data.back(); ++i)
        {
          const unsigned int dofs_mcell =
            dofs_per_cell[cell_active_fe_index.size() == 0 ? 0 :
                          cell_active_fe_index[i]] * vectorization_length;
          new_row_starts[i][0] = new_dof_indices.size();
          new_row_starts[i][1] = new_constraint_indicator.size();
          new_row_starts[i][2] = irregular_cells[i];

          const unsigned int n_comp = (irregular_cells[i]>0 ?
                                       irregular_cells[i] : vectorization_length);

          for (unsigned int j=0; j<n_comp; ++j)
            {
              glob_indices[j] = begin_indices(renumbering[position_cell+j]);
              constr_ind[j] = begin_indicators(renumbering[position_cell+j]);
              constr_end[j] = end_indicators(renumbering[position_cell+j]);
              index[j] = 0;
            }

          bool has_constraints = false;
          if (store_plain_indices == true)
            {
              for (unsigned int j=0; j<n_comp; ++j)
                if (begin_indicators(renumbering[position_cell+j]) <
                    end_indicators(renumbering[position_cell+j]))
                  {
                    plain_glob_indices[j] =
                      begin_indices_plain (renumbering[position_cell+j]);
                    has_constraints = true;
                  }
                else
                  plain_glob_indices[j] =
                    begin_indices (renumbering[position_cell+j]);
              if (has_constraints == true)
                new_rowstart_plain[i] = new_plain_indices.size();
            }

          unsigned int m_ind_local = 0, m_index = 0;
          while (m_ind_local < dofs_mcell)
            for (unsigned int j=0; j<vectorization_length; ++j)
              {
                // last cell: nothing to do
                if (j >= n_comp)
                  {
                    ++m_ind_local;
                    continue;
                  }

                // otherwise, check if we are a constrained dof. The dof is
                // not constrained if we are at the end of the row for the
                // constraints (indi[j] == n_indi[j]) or if the local index[j]
                // is smaller than the next position for a constraint. Then,
                // just copy it. otherwise, copy all the entries that come
                // with this dof
                if (constr_ind[j] == constr_end[j] ||
                    index[j] < constr_ind[j]->first)
                  {
                    new_dof_indices.push_back (*glob_indices[j]);
                    ++m_index;
                    ++index[j];
                    ++glob_indices[j];
                  }
                else
                  {
                    const unsigned short constraint_loc = constr_ind[j]->second;
                    new_constraint_indicator.emplace_back (m_index, constraint_loc);
                    for (unsigned int k=constraint_pool_row_index[constraint_loc];
                         k<constraint_pool_row_index[constraint_loc+1];
                         ++k, ++glob_indices[j])
                      new_dof_indices.push_back (*glob_indices[j]);
                    ++constr_ind[j];
                    m_index = 0;
                    index[j] = 0;
                  }
                if (store_plain_indices==true && has_constraints==true)
                  new_plain_indices.push_back (*plain_glob_indices[j]++);
                ++m_ind_local;
              }

          for (unsigned int j=0; j<n_comp; ++j)
            Assert (glob_indices[j]==end_indices(renumbering[position_cell+j]),
                    ExcInternalError());
          position_cell += n_comp;
        }
      AssertDimension (position_cell+1, row_starts.size());

      new_row_starts[task_info.cell_partition_data.back()][0] = new_dof_indices.size();
      new_row_starts[task_info.cell_partition_data.back()][1] = new_constraint_indicator.size();
      new_row_starts[task_info.cell_partition_data.back()][2] = 0;

      AssertDimension(dof_indices.size(), new_dof_indices.size());
      AssertDimension(constraint_indicator.size(),
                      new_constraint_indicator.size());

      new_row_starts.swap (row_starts);
      new_dof_indices.swap (dof_indices);
      new_constraint_indicator.swap (constraint_indicator);
      new_plain_indices.swap (plain_dof_indices);
      new_rowstart_plain.swap (row_starts_plain_indices);

#ifdef DEBUG
      // sanity check 1: all indices should be smaller than the number of dofs
      // locally owned plus the number of ghosts
      const unsigned int index_range = (vector_partitioner->local_range().second-
                                        vector_partitioner->local_range().first)
                                       + vector_partitioner->ghost_indices().n_elements();
      for (std::size_t i=0; i<dof_indices.size(); ++i)
        AssertIndexRange (dof_indices[i], index_range);

      // sanity check 2: for the constraint indicators, the first index should
      // be smaller than the number of indices in the row, and the second
      // index should be smaller than the number of constraints in the
      // constraint pool.
      for (unsigned int row=0; row<task_info.cell_partition_data.back(); ++row)
        {
          const unsigned int row_length_ind = row_length_indices(row);
          const std::pair<unsigned short,unsigned short>
          *con_it = begin_indicators(row), * end_con = end_indicators(row);
          for ( ; con_it != end_con; ++con_it)
            {
              AssertIndexRange (con_it->first, row_length_ind+1);
              AssertIndexRange (con_it->second,
                                constraint_pool_row_index.size()-1);
            }
        }

      // sanity check 3: check the number of cells once again
      unsigned int n_active_cells = 0;
      for (unsigned int c=0; c<*(task_info.cell_partition_data.end()-2); ++c)
        if (irregular_cells[c] > 0)
          n_active_cells += irregular_cells[c];
        else
          n_active_cells += vectorization_length;
      AssertDimension(n_active_cells, task_info.n_active_cells);
#endif
    }



    namespace
    {
      // rudimentary version of a vector that keeps entries always ordered
      class ordered_vector : public std::vector<types::global_dof_index>
      {
      public:
        ordered_vector ()
        {
          reserve (2000);
        }

        void reserve (const std::size_t size)
        {
          if (size > 0)
            this->std::vector<types::global_dof_index>::reserve (size);
        }


        // insert a given entry. dat is a pointer within this vector (the user
        // needs to make sure that it really stays there)
        void insert (const unsigned int entry,
                     std::vector<types::global_dof_index>::iterator &dat)
        {
          AssertIndexRange (static_cast<std::size_t>(dat - begin()), size()+1);
          AssertIndexRange (static_cast<std::size_t>(end() - dat), size()+1);
          AssertIndexRange (size(), capacity());
          while (dat != end() && *dat < entry)
            ++dat;

          if (dat == end())
            {
              push_back(entry);
              dat = end();
            }
          else if (*dat > entry)
            {
              dat = this->std::vector<types::global_dof_index>::insert (dat, entry);
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

      void compute_row_lengths(const unsigned int           begin,
                               const unsigned int           end,
                               const DoFInfo               &dof_info,
                               std::vector<Threads::Mutex> &mutexes,
                               std::vector<unsigned int>   &row_lengths)
      {
        std::vector<unsigned int> scratch;
        constexpr unsigned int n_components = 1;
        for (unsigned int block=begin; block<end; ++block)
          {
            scratch.clear();
            scratch.insert(scratch.end(),
                           &dof_info.dof_indices[dof_info.row_starts[block*n_components][0]],
                           &dof_info.dof_indices[dof_info.row_starts[(block+1)*n_components][0]]);
            std::sort(scratch.begin(), scratch.end());
            std::vector<unsigned int>::const_iterator end_unique =
              std::unique(scratch.begin(), scratch.end());
            std::vector<unsigned int>::const_iterator it = scratch.begin();
            while (it != end_unique)
              {
                // In this code, the procedure is that we insert all elements
                // that are within the range of one lock at once
                const unsigned int next_bucket = (*it/bucket_size_threading+1)*
                                                 bucket_size_threading;
                Threads::Mutex::ScopedLock lock(mutexes[*it/bucket_size_threading]);
                for ( ; it != end_unique && *it < next_bucket; ++it)
                  {
                    AssertIndexRange(*it, row_lengths.size());
                    row_lengths[*it]++;
                  }
              }
          }
      }

      void fill_connectivity_dofs(const unsigned int               begin,
                                  const unsigned int               end,
                                  const DoFInfo                   &dof_info,
                                  const std::vector<unsigned int> &row_lengths,
                                  std::vector<Threads::Mutex>     &mutexes,
                                  dealii::SparsityPattern         &connectivity_dof)
      {
        std::vector<unsigned int> scratch;
        const unsigned int n_components = 1;
        for (unsigned int block=begin; block<end; ++block)
          {
            scratch.clear();
            scratch.insert(scratch.end(),
                           &dof_info.dof_indices[dof_info.row_starts[block*n_components][0]],
                           &dof_info.dof_indices[dof_info.row_starts[(block+1)*n_components][0]]);
            std::sort(scratch.begin(), scratch.end());
            std::vector<unsigned int>::const_iterator end_unique =
              std::unique(scratch.begin(), scratch.end());
            std::vector<unsigned int>::const_iterator it = scratch.begin();
            while (it != end_unique)
              {
                const unsigned int next_bucket = (*it/bucket_size_threading+1)*
                                                 bucket_size_threading;
                Threads::Mutex::ScopedLock lock(mutexes[*it/bucket_size_threading]);
                for ( ; it != end_unique && *it < next_bucket; ++it)
                  if (row_lengths[*it]>0)
                    connectivity_dof.add(*it, block);
              }
          }
      }

      void fill_connectivity(const unsigned int               begin,
                             const unsigned int               end,
                             const DoFInfo                   &dof_info,
                             const std::vector<unsigned int> &renumbering,
                             const dealii::SparsityPattern   &connectivity_dof,
                             DynamicSparsityPattern          &connectivity)
      {
        ordered_vector row_entries;
        const unsigned int n_components = 1;
        for (unsigned int block=begin; block < end; ++block)
          {
            row_entries.clear();

            const unsigned int
            *it = &dof_info.dof_indices[dof_info.row_starts[block*n_components][0]],
             *end_cell = &dof_info.dof_indices[dof_info.row_starts[(block+1)*n_components][0]];
            for ( ; it != end_cell; ++it)
              {
                SparsityPattern::iterator sp = connectivity_dof.begin(*it);
                std::vector<types::global_dof_index>::iterator insert_pos = row_entries.begin();
                for ( ; sp != connectivity_dof.end(*it); ++sp)
                  if (sp->column() != block)
                    row_entries.insert (renumbering[sp->column()], insert_pos);
              }
            connectivity.add_entries (renumbering[block], row_entries.begin(), row_entries.end());
          }
      }
    }


    void
    DoFInfo::make_connectivity_graph
    (const TaskInfo                  &task_info,
     const std::vector<unsigned int> &renumbering,
     DynamicSparsityPattern          &connectivity) const
    {
      unsigned int n_rows =
        (vector_partitioner->local_range().second-
         vector_partitioner->local_range().first)
        + vector_partitioner->ghost_indices().n_elements();

      // Avoid square sparsity patterns that allocate the diagonal entry
      if (n_rows == task_info.n_active_cells)
        ++n_rows;

      // first determine row lengths
      std::vector<unsigned int> row_lengths(n_rows);
      std::vector<Threads::Mutex> mutexes(n_rows/bucket_size_threading+1);
      parallel::apply_to_subranges(0, task_info.n_active_cells,
                                   std::bind(&compute_row_lengths,
                                             std::placeholders::_1,
                                             std::placeholders::_2,
                                             std::cref(*this),
                                             std::ref(mutexes),
                                             std::ref(row_lengths)), 20);

      // disregard dofs that only sit on a single cell because they cannot
      // couple
      for (unsigned int row=0; row<n_rows; ++row)
        if (row_lengths[row] <= 1)
          row_lengths[row] = 0;

      // Create a temporary sparsity pattern that holds to each degree of
      // freedom on which cells it appears, i.e., store the connectivity
      // between cells and dofs
      SparsityPattern connectivity_dof (n_rows, task_info.n_active_cells,
                                        row_lengths);
      parallel::apply_to_subranges(0, task_info.n_active_cells,
                                   std::bind(&fill_connectivity_dofs,
                                             std::placeholders::_1,
                                             std::placeholders::_2,
                                             std::cref(*this),
                                             std::cref(row_lengths),
                                             std::ref(mutexes),
                                             std::ref(connectivity_dof)), 20);
      connectivity_dof.compress();


      // Invert renumbering for use in fill_connectivity.
      std::vector<unsigned int> reverse_numbering(task_info.n_active_cells);
      reverse_numbering = Utilities::invert_permutation(renumbering);

      // From the above connectivity between dofs and cells, we can finally
      // create a connectivity list between cells. The connectivity graph
      // should apply the renumbering, i.e., the entry for cell j is the entry
      // for cell renumbering[j] in the original ordering.
      parallel::apply_to_subranges(0, task_info.n_active_cells,
                                   std::bind(&fill_connectivity,
                                             std::placeholders::_1,
                                             std::placeholders::_2,
                                             std::cref(*this),
                                             std::cref(reverse_numbering),
                                             std::cref(connectivity_dof),
                                             std::ref(connectivity)), 20);
    }



    void DoFInfo::renumber_dofs (std::vector<types::global_dof_index> &renumbering)
    {
      // first renumber all locally owned degrees of freedom
      AssertDimension (vector_partitioner->local_size(),
                       vector_partitioner->size());
      const unsigned int local_size = vector_partitioner->local_size();
      renumbering.resize (0);
      renumbering.resize (local_size, numbers::invalid_dof_index);

      types::global_dof_index counter = 0;
      std::vector<unsigned int>::iterator dof_ind = dof_indices.begin(),
                                          end_ind = dof_indices.end();
      for ( ; dof_ind != end_ind; ++dof_ind)
        {
          if (*dof_ind < local_size)
            {
              if (renumbering[*dof_ind] == numbers::invalid_dof_index)
                renumbering[*dof_ind] = counter++;
              *dof_ind = renumbering[*dof_ind];
            }
        }

      AssertIndexRange (counter, local_size+1);
      for (std::size_t i=0; i<renumbering.size(); ++i)
        if (renumbering[i] == numbers::invalid_dof_index)
          renumbering[i] = counter++;

      // adjust the constrained DoFs
      std::vector<unsigned int> new_constrained_dofs (constrained_dofs.size());
      for (std::size_t i=0; i<constrained_dofs.size(); ++i)
        new_constrained_dofs[i] = renumbering[constrained_dofs[i]];

      // the new constrained DoFs should be sorted already as they are not
      // contained in dof_indices and then get contiguous numbers
#ifdef DEBUG
      for (std::size_t i=1; i<new_constrained_dofs.size(); ++i)
        Assert (new_constrained_dofs[i] > new_constrained_dofs[i-1], ExcInternalError());
#endif
      std::swap (constrained_dofs, new_constrained_dofs);

      // transform indices to global index space
      for (std::size_t i=0; i<renumbering.size(); ++i)
        renumbering[i] = vector_partitioner->local_to_global(renumbering[i]);

      AssertDimension (counter, renumbering.size());
    }



    std::size_t
    DoFInfo::memory_consumption () const
    {
      std::size_t memory = sizeof(*this);
      memory += (row_starts.capacity()*sizeof(std::array<unsigned int,3>));
      memory += MemoryConsumption::memory_consumption (dof_indices);
      memory += MemoryConsumption::memory_consumption (row_starts_plain_indices);
      memory += MemoryConsumption::memory_consumption (plain_dof_indices);
      memory += MemoryConsumption::memory_consumption (constraint_indicator);
      memory += MemoryConsumption::memory_consumption (*vector_partitioner);
      return memory;
    }



    template <typename StreamType>
    void
    DoFInfo::print_memory_consumption (StreamType     &out,
                                       const TaskInfo &task_info) const
    {
      out << "       Memory row starts indices:    ";
      task_info.print_memory_statistics
      (out, (row_starts.capacity()*sizeof(std::array<unsigned int, 3>)));
      out << "       Memory dof indices:           ";
      task_info.print_memory_statistics
      (out, MemoryConsumption::memory_consumption (dof_indices));
      out << "       Memory constraint indicators: ";
      task_info.print_memory_statistics
      (out, MemoryConsumption::memory_consumption (constraint_indicator));
      out << "       Memory plain indices:         ";
      task_info.print_memory_statistics
      (out, MemoryConsumption::memory_consumption (row_starts_plain_indices)+
       MemoryConsumption::memory_consumption (plain_dof_indices));
      out << "       Memory vector partitioner:    ";
      task_info.print_memory_statistics
      (out, MemoryConsumption::memory_consumption (*vector_partitioner));
    }



    template <typename Number>
    void
    DoFInfo::print (const std::vector<Number>       &constraint_pool_data,
                    const std::vector<unsigned int> &constraint_pool_row_index,
                    std::ostream                    &out) const
    {
      const unsigned int n_rows = row_starts.size() - 1;
      for (unsigned int row=0 ; row<n_rows ; ++row)
        {
          out << "Entries row " << row << ": ";
          const unsigned int *glob_indices = begin_indices(row),
                              *end_row = end_indices(row);
          unsigned int index = 0;
          const std::pair<unsigned short,unsigned short>
          *con_it = begin_indicators(row),
           * end_con = end_indicators(row);
          for ( ; con_it != end_con; ++con_it)
            {
              for ( ; index<con_it->first; index++)
                {
                  Assert (glob_indices+index != end_row, ExcInternalError());
                  out << glob_indices[index] << " ";
                }

              out << "[ ";
              for (unsigned int k=constraint_pool_row_index[con_it->second];
                   k<constraint_pool_row_index[con_it->second+1];
                   k++,index++)
                {
                  Assert (glob_indices+index != end_row, ExcInternalError());
                  out << glob_indices[index] << "/"
                      << constraint_pool_data[k];
                  if (k<constraint_pool_row_index[con_it->second+1]-1)
                    out << " ";
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
