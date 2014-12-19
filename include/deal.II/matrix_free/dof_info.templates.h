// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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


#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/helper_functions.h>

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
      // equal. this was quite lenghty and now we use a std::map with a
      // user-defined comparator to compare floating point arrays to a
      // tolerance 1e-13.
      std::pair<typename std::map<std::vector<double>, types::global_dof_index,
          FPArrayComparator<double> >::iterator,
          bool> it = constraints.insert(next_constraint);

      types::global_dof_index insert_position = deal_II_numbers::invalid_dof_index;
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


    DoFInfo::DoFInfo (const DoFInfo &dof_info_in)
      :
      row_starts (dof_info_in.row_starts),
      dof_indices (dof_info_in.dof_indices),
      constraint_indicator (dof_info_in.constraint_indicator),
      vector_partitioner (dof_info_in.vector_partitioner),
      constrained_dofs (dof_info_in.constrained_dofs),
      row_starts_plain_indices (dof_info_in.row_starts_plain_indices),
      plain_dof_indices (dof_info_in.plain_dof_indices),
      dimension (dof_info_in.dimension),
      n_components (dof_info_in.n_components),
      dofs_per_cell (dof_info_in.dofs_per_cell),
      dofs_per_face (dof_info_in.dofs_per_face),
      store_plain_indices (dof_info_in.store_plain_indices),
      cell_active_fe_index (dof_info_in.cell_active_fe_index),
      max_fe_index (dof_info_in.max_fe_index),
      fe_index_conversion (dof_info_in.fe_index_conversion),
      ghost_dofs (dof_info_in.ghost_dofs)
    {}



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
      Assert (vector_partitioner.get() !=0, ExcInternalError());
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
          if (entries_ptr != 0)
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
      unsigned int      n_unique_ghosts= 0;
#ifdef DEBUG
      for (std::vector<unsigned int>::iterator dof = dof_indices.begin();
           dof!=dof_indices.end(); ++dof)
        AssertIndexRange (*dof, n_owned+n_ghosts);
#endif

      std::vector<unsigned int> ghost_numbering (n_ghosts);
      IndexSet ghost_indices (vector_partitioner->size());
      if (n_ghosts > 0)
        {
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
    DoFInfo::compute_renumber_serial (const std::vector<unsigned int> &boundary_cells,
                                      const SizeInfo                  &size_info,
                                      std::vector<unsigned int>       &renumbering)
    {
      std::vector<unsigned int> reverse_numbering (size_info.n_active_cells,
                                                   numbers::invalid_unsigned_int);
      const unsigned int n_boundary_cells = boundary_cells.size();
      for (unsigned int j=0; j<n_boundary_cells; ++j)
        reverse_numbering[boundary_cells[j]] =
          j + size_info.vectorization_length*size_info.boundary_cells_start;
      unsigned int counter = 0;
      unsigned int j = 0;
      while (counter < size_info.n_active_cells &&
             counter < size_info.vectorization_length * size_info.boundary_cells_start)
        {
          if (reverse_numbering[j] == numbers::invalid_unsigned_int)
            reverse_numbering[j] = counter++;
          j++;
        }
      counter = std::min (size_info.vectorization_length*
                          size_info.boundary_cells_start+n_boundary_cells,
                          size_info.n_active_cells);
      if (counter < size_info.n_active_cells)
        {
          for ( ; j<size_info.n_active_cells; ++j)
            if (reverse_numbering[j] == numbers::invalid_unsigned_int)
              reverse_numbering[j] = counter++;
        }
      AssertDimension (counter, size_info.n_active_cells);
      renumbering = Utilities::invert_permutation (reverse_numbering);
    }



    void
    DoFInfo::compute_renumber_hp_serial (SizeInfo                  &size_info,
                                         std::vector<unsigned int> &renumbering,
                                         std::vector<unsigned int> &irregular_cells)
    {
      if (max_fe_index < 2)
        return;
      const unsigned int n_active_cells = size_info.n_active_cells;
      const unsigned int vectorization_length = size_info.vectorization_length;
      irregular_cells.resize (0);
      irregular_cells.resize (size_info.n_macro_cells+3*max_fe_index);
      std::vector<std::vector<unsigned int> > renumbering_fe_index;
      renumbering_fe_index.resize(max_fe_index);
      unsigned int counter,n_macro_cells_before = 0;
      const unsigned int
      start_bound = std::min (size_info.n_active_cells,
                              size_info.boundary_cells_start*vectorization_length),
                    end_bound   = std::min (size_info.n_active_cells,
                                            size_info.boundary_cells_end*vectorization_length);
      for (counter=0; counter<start_bound; counter++)
        {
          renumbering_fe_index[cell_active_fe_index[renumbering[counter]]].
          push_back(renumbering[counter]);
        }
      counter = 0;
      for (unsigned int j=0; j<max_fe_index; j++)
        {
          for (unsigned int jj=0; jj<renumbering_fe_index[j].size(); jj++)
            renumbering[counter++] = renumbering_fe_index[j][jj];
          irregular_cells[renumbering_fe_index[j].size()/vectorization_length+
                          n_macro_cells_before] =
                            renumbering_fe_index[j].size()%vectorization_length;
          n_macro_cells_before += (renumbering_fe_index[j].size()+vectorization_length-1)/
                                  vectorization_length;
          renumbering_fe_index[j].resize(0);
        }
      unsigned int new_boundary_start = n_macro_cells_before;
      for (counter = start_bound; counter < end_bound; counter++)
        {
          renumbering_fe_index[cell_active_fe_index[renumbering[counter]]].
          push_back(renumbering[counter]);
        }
      counter = start_bound;
      for (unsigned int j=0; j<max_fe_index; j++)
        {
          for (unsigned int jj=0; jj<renumbering_fe_index[j].size(); jj++)
            renumbering[counter++] = renumbering_fe_index[j][jj];
          irregular_cells[renumbering_fe_index[j].size()/vectorization_length+
                          n_macro_cells_before] =
                            renumbering_fe_index[j].size()%vectorization_length;
          n_macro_cells_before += (renumbering_fe_index[j].size()+vectorization_length-1)/
                                  vectorization_length;
          renumbering_fe_index[j].resize(0);
        }
      unsigned int new_boundary_end = n_macro_cells_before;
      for (counter=end_bound; counter<n_active_cells; counter++)
        {
          renumbering_fe_index[cell_active_fe_index[renumbering[counter]]].
          push_back(renumbering[counter]);
        }
      counter = end_bound;
      for (unsigned int j=0; j<max_fe_index; j++)
        {
          for (unsigned int jj=0; jj<renumbering_fe_index[j].size(); jj++)
            renumbering[counter++] = renumbering_fe_index[j][jj];
          irregular_cells[renumbering_fe_index[j].size()/vectorization_length+
                          n_macro_cells_before] =
                            renumbering_fe_index[j].size()%vectorization_length;
          n_macro_cells_before += (renumbering_fe_index[j].size()+vectorization_length-1)/
                                  vectorization_length;
        }
      AssertIndexRange (n_macro_cells_before,
                        size_info.n_macro_cells + 3*max_fe_index+1);
      irregular_cells.resize (n_macro_cells_before);
      size_info.n_macro_cells = n_macro_cells_before;
      size_info.boundary_cells_start = new_boundary_start;
      size_info.boundary_cells_end = new_boundary_end;
    }



    void
    DoFInfo::compute_renumber_parallel (const std::vector<unsigned int> &boundary_cells,
                                        SizeInfo                        &size_info,
                                        std::vector<unsigned int>       &renumbering)
    {
      std::vector<unsigned int> reverse_numbering (size_info.n_active_cells,
                                                   numbers::invalid_unsigned_int);
      const unsigned int n_boundary_cells = boundary_cells.size();
      for (unsigned int j=0; j<n_boundary_cells; ++j)
        reverse_numbering[boundary_cells[j]] = j;
      unsigned int counter = n_boundary_cells;
      for (unsigned int j=0; j<size_info.n_active_cells; ++j)
        if (reverse_numbering[j] == numbers::invalid_unsigned_int)
          reverse_numbering[j] = counter++;

      size_info.boundary_cells_end   = (size_info.boundary_cells_end -
                                        size_info.boundary_cells_start);
      size_info.boundary_cells_start = 0;

      AssertDimension (counter, size_info.n_active_cells);
      renumbering = Utilities::invert_permutation (reverse_numbering);
    }



    void
    DoFInfo::reorder_cells (const SizeInfo                  &size_info,
                            const std::vector<unsigned int> &renumbering,
                            const std::vector<unsigned int> &constraint_pool_row_index,
                            const std::vector<unsigned int> &irregular_cells,
                            const unsigned int               vectorization_length)
    {
      // first reorder the active fe index.
      if (cell_active_fe_index.size() > 0)
        {
          std::vector<unsigned int> new_active_fe_index;
          new_active_fe_index.reserve (size_info.n_macro_cells);
          std::vector<unsigned int> fe_indices(vectorization_length);
          unsigned int position_cell = 0;
          for (unsigned int cell=0; cell<size_info.n_macro_cells; ++cell)
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

      std::vector<std_cxx11::array<unsigned int, 3> > new_row_starts;
      std::vector<unsigned int> new_dof_indices;
      std::vector<std::pair<unsigned short,unsigned short> >
      new_constraint_indicator;
      std::vector<unsigned int> new_plain_indices, new_rowstart_plain;
      unsigned int position_cell = 0;
      new_row_starts.resize (size_info.n_macro_cells + 1);
      new_dof_indices.reserve (dof_indices.size());
      new_constraint_indicator.reserve (constraint_indicator.size());
      if (store_plain_indices == true)
        {
          new_rowstart_plain.resize (size_info.n_macro_cells + 1,
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
      for (unsigned int i=0; i<size_info.n_macro_cells; ++i)
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
                    new_constraint_indicator.push_back
                    (std::pair<unsigned short,unsigned short> (m_index, constraint_loc));
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

      new_row_starts[size_info.n_macro_cells][0] = new_dof_indices.size();
      new_row_starts[size_info.n_macro_cells][1] = new_constraint_indicator.size();
      new_row_starts[size_info.n_macro_cells][2] = 0;

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
      for (unsigned int row=0; row<size_info.n_macro_cells; ++row)
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

      // sanity check 3: all non-boundary cells should have indices that only
      // refer to the locally owned range
      const unsigned int local_size = (vector_partitioner->local_range().second-
                                       vector_partitioner->local_range().first);
      for (unsigned int row=0; row<size_info.boundary_cells_start; ++row)
        {
          const unsigned int *ptr     = begin_indices(row);
          const unsigned int *end_ptr = end_indices  (row);
          for ( ; ptr != end_ptr; ++ptr)
            AssertIndexRange (*ptr, local_size);
        }
      for (unsigned int row=size_info.boundary_cells_end;
           row<size_info.n_macro_cells; ++row)
        {
          const unsigned int *ptr     = begin_indices(row);
          const unsigned int *end_ptr = end_indices  (row);
          for ( ; ptr != end_ptr; ++ptr)
            AssertIndexRange (*ptr, local_size);
        }
#endif
    }



    void DoFInfo::guess_block_size (const SizeInfo &size_info,
                                    TaskInfo       &task_info)
    {
      // user did not say a positive number, so we have to guess
      if (task_info.block_size == 0)
        {
          // we would like to have enough work to do, so as first guess, try
          // to get 50 times as many chunks as we have threads on the system.
          task_info.block_size =
            size_info.n_macro_cells / (multithread_info.n_threads() * 50);

          // if there are too few degrees of freedom per cell, need to
          // increase the block size
          const unsigned int minimum_parallel_grain_size = 500;
          if (dofs_per_cell[0] * task_info.block_size <
              minimum_parallel_grain_size)
            task_info.block_size = (minimum_parallel_grain_size /
                                    dofs_per_cell[0] + 1);
        }
      if (task_info.block_size > size_info.n_macro_cells)
        task_info.block_size = size_info.n_macro_cells;
    }



    void DoFInfo::make_thread_graph_partition_color
    (SizeInfo                  &size_info,
     TaskInfo                  &task_info,
     std::vector<unsigned int> &renumbering,
     std::vector<unsigned int> &irregular_cells,
     const bool                 hp_bool)
    {
      if (size_info.n_macro_cells == 0)
        return;

      const std::size_t vectorization_length = size_info.vectorization_length;
      Assert (vectorization_length > 0, ExcInternalError());

      guess_block_size (size_info, task_info);

      // set up partitions. if we just use coloring without partitions, do
      // nothing here, assume all cells to belong to the zero partition (that
      // we otherwise use for MPI boundary cells)
      unsigned int start_up = 0,
                   start_nonboundary = numbers::invalid_unsigned_int;
      if (task_info.use_coloring_only == false)
        {
          start_nonboundary =
            std::min(((size_info.boundary_cells_end+task_info.block_size-1)/
                      task_info.block_size)*task_info.block_size,
                     size_info.n_macro_cells);
          start_up = start_nonboundary;
          size_info.boundary_cells_end = start_nonboundary;
        }
      else
        {
          start_nonboundary = size_info.n_macro_cells;
          start_up = size_info.n_macro_cells;
          size_info.boundary_cells_start = 0;
          size_info.boundary_cells_end = size_info.n_macro_cells;
        }
      if (hp_bool == true)
        {
          irregular_cells.resize (0);
          irregular_cells.resize (size_info.n_macro_cells+2*max_fe_index);
          std::vector<std::vector<unsigned int> > renumbering_fe_index;
          renumbering_fe_index.resize(max_fe_index);
          unsigned int counter,n_macro_cells_before = 0;
          for (counter=0; counter<start_nonboundary*vectorization_length;
               counter++)
            {
              renumbering_fe_index[cell_active_fe_index[renumbering[counter]]].
              push_back(renumbering[counter]);
            }
          counter = 0;
          for (unsigned int j=0; j<max_fe_index; j++)
            {
              for (unsigned int jj=0; jj<renumbering_fe_index[j].size(); jj++)
                renumbering[counter++] = renumbering_fe_index[j][jj];
              irregular_cells[renumbering_fe_index[j].size()/vectorization_length+
                              n_macro_cells_before] =
                                renumbering_fe_index[j].size()%vectorization_length;
              n_macro_cells_before += (renumbering_fe_index[j].size()+vectorization_length-1)/
                                      vectorization_length;
              renumbering_fe_index[j].resize(0);
            }

          unsigned int new_boundary_end = n_macro_cells_before;
          for (counter=start_nonboundary*vectorization_length;
               counter<size_info.n_active_cells; counter++)
            {
              renumbering_fe_index[cell_active_fe_index.empty() ? 0 :
                                   cell_active_fe_index[renumbering[counter]]].
              push_back(renumbering[counter]);
            }
          counter = start_nonboundary * vectorization_length;
          for (unsigned int j=0; j<max_fe_index; j++)
            {
              for (unsigned int jj=0; jj<renumbering_fe_index[j].size(); jj++)
                renumbering[counter++] = renumbering_fe_index[j][jj];
              irregular_cells[renumbering_fe_index[j].size()/vectorization_length+
                              n_macro_cells_before] =
                                renumbering_fe_index[j].size()%vectorization_length;
              n_macro_cells_before += (renumbering_fe_index[j].size()+vectorization_length-1)/
                                      vectorization_length;
            }
          AssertIndexRange (n_macro_cells_before,
                            size_info.n_macro_cells + 2*max_fe_index+1);
          irregular_cells.resize (n_macro_cells_before);
          size_info.n_macro_cells = n_macro_cells_before;
          size_info.boundary_cells_start = 0;
          size_info.boundary_cells_end = new_boundary_end;
          task_info.n_blocks = (size_info.n_macro_cells+task_info.block_size-1)
                               /task_info.block_size;
          task_info.block_size_last = size_info.n_macro_cells%task_info.block_size;
          if (task_info.block_size_last == 0)
            task_info.block_size_last = task_info.block_size;
        }

      // assume that all FEs have the same connectivity graph, so take the
      // zeroth FE
      task_info.n_blocks = (size_info.n_macro_cells+task_info.block_size-1)/
                           task_info.block_size;
      task_info.block_size_last = size_info.n_macro_cells-
                                  (task_info.block_size*(task_info.n_blocks-1));

      // create the connectivity graph with internal blocking
      CompressedSimpleSparsityPattern connectivity;
      make_connectivity_graph (size_info, task_info, renumbering,irregular_cells,
                               true, connectivity);

      // Create cell-block  partitioning.
      unsigned int partition = 0, counter = 0;
      bool work = true;

      // For each block of cells, this variable saves to which partitions the
      // block belongs. Initialize all to n_macro_cells to mark them as not
      // yet assigned a partition.
      std::vector<unsigned int> cell_partition(task_info.n_blocks,
                                               size_info.n_macro_cells);
      std::vector<unsigned int> neighbor_list;
      std::vector<unsigned int> neighbor_neighbor_list;

      // In element j of this variable, one puts the old number of the block
      // that should be the jth block in the new numeration.
      std::vector<unsigned int> partition_list      (task_info.n_blocks,0);
      std::vector<unsigned int> partition_color_list(task_info.n_blocks,0);

      // This vector points to the start of each partition.
      std::vector<unsigned int> partition_blocks (2,0);
      std::vector<unsigned int> cell_color(task_info.n_blocks,
                                           size_info.n_macro_cells);
      std::vector<bool> color_finder;

      // this performs a classical breath-first search in the connectivity
      // graph of the cell chunks
      while (work)
        {
          // put all cells up to begin_inner_cells into first partition. if
          // the numbers do not add up exactly, assign an additional block
          if (start_nonboundary>0 && start_up == start_nonboundary)
            {
              unsigned int n_blocks = ((start_nonboundary+task_info.block_size-1)
                                       /task_info.block_size);
              for (unsigned int cell=0; cell<n_blocks; ++cell)
                {
                  cell_partition[cell] = partition;
                  neighbor_list.push_back(cell);
                  partition_list[counter++] = cell;
                  partition_blocks.back()++;
                }
            }
          else
            {
              // To start up, set the start_up cell to partition and list all
              // its neighbors.
              AssertIndexRange(start_up, cell_partition.size());
              cell_partition[start_up] = partition;
              neighbor_list.push_back(start_up);
              partition_list[counter++] = start_up;
              partition_blocks.back()++;
            }

          while (neighbor_list.size()>0)
            {
              partition++;
              partition_blocks.push_back(partition_blocks.back());
              for (unsigned int j=0; j<neighbor_list.size(); ++j)
                {
                  Assert(cell_partition[neighbor_list[j]]==partition-1,
                         ExcInternalError());
                  CompressedSimpleSparsityPattern::row_iterator neighbor =
                    connectivity.row_begin(neighbor_list[j]),
                    end = connectivity.row_end(neighbor_list[j]);
                  for (; neighbor!=end ; ++neighbor)
                    {
                      if (cell_partition[*neighbor]==size_info.n_macro_cells)
                        {
                          partition_blocks.back()++;
                          cell_partition[*neighbor] = partition;
                          neighbor_neighbor_list.push_back(*neighbor);
                          partition_list[counter++] = *neighbor;
                        }
                    }
                }
              neighbor_list = neighbor_neighbor_list;
              neighbor_neighbor_list.resize(0);
            }

          // One has to check if the graph is not connected so we have to find
          // another partition.
          work = false;
          for (unsigned int j=start_up; j<task_info.n_blocks; ++j)
            if (cell_partition[j] == size_info.n_macro_cells)
              {
                start_up = j;
                work = true;
                break;
              }
        }
      AssertDimension (partition_blocks[partition], task_info.n_blocks);


      // Color the cells within each partition
      task_info.partition_color_blocks_row_index.resize(partition+1);
      unsigned int color_counter = 0, index_counter = 0;
      for (unsigned int part=0; part<partition; part++)
        {
          task_info.partition_color_blocks_row_index[part] = index_counter;
          unsigned int max_color = 0;
          for (unsigned int k=partition_blocks[part]; k<partition_blocks[part+1];
               k++)
            {
              unsigned int cell = partition_list[k];
              unsigned int n_neighbors = connectivity.row_length(cell);

              // In the worst case, each neighbor has a different color. So we
              // find at least one available color between 0 and n_neighbors.
              color_finder.resize(n_neighbors+1);
              for (unsigned int j=0; j<=n_neighbors; ++j)
                color_finder[j]=true;
              CompressedSimpleSparsityPattern::row_iterator
              neighbor = connectivity.row_begin(cell),
              end      = connectivity.row_end(cell);
              for (; neighbor!=end ; ++neighbor)
                {
                  // Mark the color that a neighbor within the partition has
                  // as taken
                  if (cell_partition[*neighbor] == part &&
                      cell_color[*neighbor] <= n_neighbors)
                    color_finder[cell_color[*neighbor]] = false;
                }
              // Choose the smallest color that is not taken for the block
              cell_color[cell]=0;
              while (color_finder[cell_color[cell]] == false)
                cell_color[cell]++;
              if (cell_color[cell] > max_color)
                max_color = cell_color[cell];
            }
          // Reorder within partition: First, all blocks that belong the 0 and
          // then so on until those with color max (Note that the smaller the
          // number the larger the partition)
          for (unsigned int color=0; color<=max_color; color++)
            {
              task_info.partition_color_blocks_data.push_back(color_counter);
              index_counter++;
              for (unsigned int k=partition_blocks[part];
                   k<partition_blocks[part+1]; k++)
                {
                  unsigned int cell=partition_list[k];
                  if (cell_color[cell] == color)
                    {
                      partition_color_list[color_counter++] = cell;
                    }
                }
            }
        }
      task_info.partition_color_blocks_data.push_back(task_info.n_blocks);
      task_info.partition_color_blocks_row_index[partition] = index_counter;
      AssertDimension (color_counter, task_info.n_blocks);

      partition_list = renumbering;

      // in debug mode, check that the partition color list is one-to-one
#ifdef DEBUG
      {
        std::vector<unsigned int> sorted_pc_list (partition_color_list);
        std::sort(sorted_pc_list.begin(), sorted_pc_list.end());
        for (unsigned int i=0; i<sorted_pc_list.size(); ++i)
          Assert(sorted_pc_list[i] == i, ExcInternalError());
      }
#endif

      // set the start list for each block and compute the renumbering of
      // cells
      std::vector<unsigned int> block_start(size_info.n_macro_cells+1);
      std::vector<unsigned int> irregular(size_info.n_macro_cells);

      unsigned int mcell_start=0;
      block_start[0] = 0;
      for (unsigned int block=0; block<task_info.n_blocks; block++)
        {
          block_start[block+1] = block_start[block];
          for (unsigned int mcell=mcell_start; mcell<
               std::min(mcell_start+task_info.block_size,
                        size_info.n_macro_cells);
               ++mcell)
            {
              unsigned int n_comp = (irregular_cells[mcell]>0)
                                    ?irregular_cells[mcell]:size_info.vectorization_length;
              block_start[block+1] += n_comp;
              ++counter;
            }
          mcell_start += task_info.block_size;
        }
      counter = 0;
      unsigned int counter_macro = 0;
      for (unsigned int block=0; block<task_info.n_blocks; block++)
        {
          unsigned int present_block = partition_color_list[block];
          for (unsigned int cell = block_start[present_block];
               cell<block_start[present_block+1]; ++cell)
            renumbering[counter++] = partition_list[cell];
          unsigned int this_block_size = (present_block == task_info.n_blocks-1)?
                                         task_info.block_size_last:task_info.block_size;
          for (unsigned int j=0; j<this_block_size; j++)
            irregular[counter_macro++] =
              irregular_cells[present_block*task_info.block_size+j];
          if (present_block == task_info.n_blocks-1)
            task_info.position_short_block = block;
        }
      irregular_cells.swap(irregular);
      AssertDimension (counter, size_info.n_active_cells);
      AssertDimension (counter_macro, size_info.n_macro_cells);

      // check that the renumbering is one-to-one
#ifdef DEBUG
      {
        std::vector<unsigned int> sorted_renumbering (renumbering);
        std::sort(sorted_renumbering.begin(), sorted_renumbering.end());
        for (unsigned int i=0; i<sorted_renumbering.size(); ++i)
          Assert(sorted_renumbering[i] == i, ExcInternalError());
      }
#endif
      AssertDimension(counter,size_info.n_active_cells);
      task_info.evens = (partition+1)/2;
      task_info.odds  = (partition)/2;
      task_info.n_blocked_workers = task_info.odds-
                                    (task_info.odds+task_info.evens+1)%2;
      task_info.n_workers = task_info.partition_color_blocks_data.size()-1-
                            task_info.n_blocked_workers;
    }



    void
    DoFInfo::make_thread_graph_partition_partition
    (SizeInfo                  &size_info,
     TaskInfo                  &task_info,
     std::vector<unsigned int> &renumbering,
     std::vector<unsigned int> &irregular_cells,
     const bool                 hp_bool)
    {
      if (size_info.n_macro_cells == 0)
        return;

      const std::size_t vectorization_length = size_info.vectorization_length;
      Assert (vectorization_length > 0, ExcInternalError());

      guess_block_size (size_info, task_info);

      // assume that all FEs have the same connectivity graph, so take the
      // zeroth FE
      task_info.n_blocks = (size_info.n_macro_cells+task_info.block_size-1)/
                           task_info.block_size;
      task_info.block_size_last = size_info.n_macro_cells-
                                  (task_info.block_size*(task_info.n_blocks-1));
      task_info.position_short_block = task_info.n_blocks-1;
      unsigned int cluster_size = task_info.block_size*vectorization_length;

      // create the connectivity graph without internal blocking
      CompressedSimpleSparsityPattern connectivity;
      make_connectivity_graph (size_info, task_info, renumbering,irregular_cells,
                               false, connectivity);

      // Create cell-block  partitioning.

      // For each block of cells, this variable saves to which partitions the
      // block belongs. Initialize all to n_macro_cells to mark them as not
      // yet assigned a partition.
      std::vector<unsigned int> cell_partition (size_info.n_active_cells,
                                                size_info.n_active_cells);
      std::vector<unsigned int> neighbor_list;
      std::vector<unsigned int> neighbor_neighbor_list;

      // In element j of this variable, one puts the old number of the block
      // that should be the jth block in the new numeration.
      std::vector<unsigned int> partition_list(size_info.n_active_cells,0);
      std::vector<unsigned int> partition_partition_list(size_info.n_active_cells,0);

      // This vector points to the start of each partition.
      std::vector<unsigned int> partition_size(2,0);

      unsigned int partition = 0,start_up=0,counter=0;
      unsigned int start_nonboundary = vectorization_length * size_info.boundary_cells_end;
      if (start_nonboundary > size_info.n_active_cells)
        start_nonboundary = size_info.n_active_cells;
      bool work = true;
      unsigned int remainder = cluster_size;

      // this performs a classical breath-first search in the connectivity
      // graph of the cells under the restriction that the size of the
      // partitions should be a multiple of the given block size
      while (work)
        {
          // put the cells with neighbors on remote MPI processes up front
          if (start_nonboundary>0)
            {
              for (unsigned int cell=0; cell<start_nonboundary; ++cell)
                {
                  const unsigned int cell_nn = renumbering[cell];
                  cell_partition[cell_nn] = partition;
                  neighbor_list.push_back(cell_nn);
                  partition_list[counter++] = cell_nn;
                  partition_size.back()++;
                }
              remainder -= (start_nonboundary%cluster_size);
              if (remainder == cluster_size)
                remainder = 0;

              // adjust end of boundary cells to the remainder
              size_info.boundary_cells_end += (remainder+vectorization_length-1)/vectorization_length;
            }
          else
            {
              // To start up, set the start_up cell to partition and list all
              // its neighbors.
              cell_partition[start_up] = partition;
              neighbor_list.push_back(start_up);
              partition_list[counter++] = start_up;
              partition_size.back()++;
              start_up++;
              remainder--;
              if (remainder == cluster_size)
                remainder = 0;
            }
          int index_before = neighbor_list.size(), index = index_before,
              index_stop = 0;
          while (remainder>0)
            {
              if (index==index_stop)
                {
                  index = neighbor_list.size();
                  if (index == index_before)
                    {
                      neighbor_list.resize(0);
                      goto not_connect;
                    }
                  index_stop = index_before;
                  index_before = index;
                }
              index--;
              unsigned int additional = neighbor_list[index];
              CompressedSimpleSparsityPattern::row_iterator neighbor =
                connectivity.row_begin(additional),
                end = connectivity.row_end(additional);
              for (; neighbor!=end ; ++neighbor)
                {
                  if (cell_partition[*neighbor]==size_info.n_active_cells)
                    {
                      partition_size.back()++;
                      cell_partition[*neighbor] = partition;
                      neighbor_list.push_back(*neighbor);
                      partition_list[counter++] = *neighbor;
                      remainder--;
                      if (remainder == 0)
                        break;
                    }
                }
            }

          while (neighbor_list.size()>0)
            {
              partition++;
              unsigned int partition_counter = 0;
              partition_size.push_back(partition_size.back());

              for (unsigned int j=0; j<neighbor_list.size(); ++j)
                {
                  Assert(cell_partition[neighbor_list[j]]==partition-1,
                         ExcInternalError());
                  CompressedSimpleSparsityPattern::row_iterator neighbor =
                    connectivity.row_begin(neighbor_list[j]),
                    end = connectivity.row_end(neighbor_list[j]);
                  for (; neighbor!=end ; ++neighbor)
                    {
                      if (cell_partition[*neighbor]==size_info.n_active_cells)
                        {
                          partition_size.back()++;
                          cell_partition[*neighbor] = partition;
                          neighbor_neighbor_list.push_back(*neighbor);
                          partition_list[counter++] = *neighbor;
                          partition_counter++;
                        }
                    }
                }
              remainder = cluster_size-(partition_counter%cluster_size);
              if (remainder == cluster_size)
                remainder = 0;
              int index_stop = 0;
              int index_before = neighbor_neighbor_list.size(), index = index_before;
              while (remainder>0)
                {
                  if (index==index_stop)
                    {
                      index = neighbor_neighbor_list.size();
                      if (index == index_before)
                        {
                          neighbor_neighbor_list.resize(0);
                          break;
                        }
                      index_stop = index_before;
                      index_before = index;
                    }
                  index--;
                  unsigned int additional = neighbor_neighbor_list[index];
                  CompressedSimpleSparsityPattern::row_iterator neighbor =
                    connectivity.row_begin(additional),
                    end = connectivity.row_end(additional);
                  for (; neighbor!=end ; ++neighbor)
                    {
                      if (cell_partition[*neighbor]==size_info.n_active_cells)
                        {
                          partition_size.back()++;
                          cell_partition[*neighbor] = partition;
                          neighbor_neighbor_list.push_back(*neighbor);
                          partition_list[counter++] = *neighbor;
                          remainder--;
                          if (remainder == 0)
                            break;
                        }
                    }
                }

              neighbor_list = neighbor_neighbor_list;
              neighbor_neighbor_list.resize(0);
            }
not_connect:
          // One has to check if the graph is not connected so we have to find
          // another partition.
          work = false;
          for (unsigned int j=start_up; j<size_info.n_active_cells; ++j)
            if (cell_partition[j] == size_info.n_active_cells)
              {
                start_up = j;
                work = true;
                if (remainder == 0)
                  remainder = cluster_size;
                break;
              }
        }
      if (remainder != 0)
        partition++;

      for (unsigned int j=0; j<renumbering.size(); j++)
        renumbering[j] = 0;
      irregular_cells.back() = 0;
      irregular_cells.resize(size_info.n_active_cells);
      unsigned int n_macro_cells_before = 0;
      {
        // Create partitioning within partitions.

        // For each block of cells, this variable saves to which partitions
        // the block belongs. Initialize all to n_macro_cells to mark them as
        // not yet assigned a partition.
        std::vector<unsigned int> cell_partition_l2(size_info.n_active_cells,
                                                    size_info.n_active_cells);
        task_info.partition_color_blocks_row_index.resize(partition+1,0);
        task_info.partition_color_blocks_data.resize(1,0);

        start_up = 0;
        counter = 0;
        unsigned int missing_macros;
        for (unsigned int part=0; part<partition; ++part)
          {
            neighbor_neighbor_list.resize(0);
            neighbor_list.resize(0);
            bool work = true;
            unsigned int partition_l2 = 0;
            start_up = partition_size[part];
            unsigned int partition_counter = 0;
            while (work)
              {
                if (neighbor_list.size()==0)
                  {
                    work = false;
                    partition_counter = 0;
                    for (unsigned int j=start_up; j<partition_size[part+1]; ++j)
                      if (cell_partition[partition_list[j]] == part &&
                          cell_partition_l2[partition_list[j]] == size_info.n_active_cells)
                        {
                          start_up = j;
                          work = true;
                          partition_counter = 1;
                          // To start up, set the start_up cell to partition
                          // and list all its neighbors.
                          AssertIndexRange (start_up, partition_size[part+1]);
                          cell_partition_l2[partition_list[start_up]] =
                            partition_l2;
                          neighbor_neighbor_list.push_back
                          (partition_list[start_up]);
                          partition_partition_list[counter++] =
                            partition_list[start_up];
                          start_up++;
                          break;
                        }
                  }
                else
                  {
                    partition_counter = 0;
                    for (unsigned int j=0; j<neighbor_list.size(); ++j)
                      {
                        Assert(cell_partition[neighbor_list[j]]==part,
                               ExcInternalError());
                        Assert(cell_partition_l2[neighbor_list[j]]==partition_l2-1,
                               ExcInternalError());
                        CompressedSimpleSparsityPattern::row_iterator neighbor =
                          connectivity.row_begin(neighbor_list[j]),
                          end = connectivity.row_end(neighbor_list[j]);
                        for (; neighbor!=end ; ++neighbor)
                          {
                            if (cell_partition[*neighbor] == part &&
                                cell_partition_l2[*neighbor]==
                                size_info.n_active_cells)
                              {
                                cell_partition_l2[*neighbor] = partition_l2;
                                neighbor_neighbor_list.push_back(*neighbor);
                                partition_partition_list[counter++] = *neighbor;
                                partition_counter++;
                              }
                          }
                      }
                  }
                if (partition_counter>0)
                  {
                    int index_before = neighbor_neighbor_list.size(),
                        index = index_before;
                    {
                      // put the cells into separate lists for each FE index
                      // within one partition-partition
                      missing_macros = 0;
                      std::vector<unsigned int> remaining_per_macro_cell
                      (max_fe_index);
                      std::vector<std::vector<unsigned int> >
                      renumbering_fe_index;
                      unsigned int cell;
                      bool filled = true;
                      if (hp_bool == true)
                        {
                          renumbering_fe_index.resize(max_fe_index);
                          for (cell=counter-partition_counter; cell<counter; ++cell)
                            {
                              renumbering_fe_index
                              [cell_active_fe_index.empty() ? 0 :
                               cell_active_fe_index[partition_partition_list
                                                    [cell]]].
                              push_back(partition_partition_list[cell]);
                            }
                          // check how many more cells are needed in the lists
                          for (unsigned int j=0; j<max_fe_index; j++)
                            {
                              remaining_per_macro_cell[j] =
                                renumbering_fe_index[j].size()%vectorization_length;
                              if (remaining_per_macro_cell[j] != 0)
                                filled = false;
                              missing_macros += ((renumbering_fe_index[j].size()+
                                                  vectorization_length-1)/vectorization_length);
                            }
                        }
                      else
                        {
                          remaining_per_macro_cell.resize(1);
                          remaining_per_macro_cell[0] = partition_counter%
                                                        vectorization_length;
                          missing_macros = partition_counter/vectorization_length;
                          if (remaining_per_macro_cell[0] != 0)
                            {
                              filled = false;
                              missing_macros++;
                            }
                        }
                      missing_macros = task_info.block_size -
                                       (missing_macros%task_info.block_size);

                      // now we realized that there are some cells missing.
                      while (missing_macros>0 || filled == false)
                        {
                          if (index==0)
                            {
                              index = neighbor_neighbor_list.size();
                              if (index == index_before)
                                {
                                  if (missing_macros != 0)
                                    {
                                      neighbor_neighbor_list.resize(0);
                                    }
                                  start_up--;
                                  break;// not connected - start again
                                }
                              index_before = index;
                            }
                          index--;
                          unsigned int additional = neighbor_neighbor_list
                                                    [index];

                          // go through the neighbors of the last cell in the
                          // current partition and check if we find some to
                          // fill up with.
                          CompressedSimpleSparsityPattern::row_iterator
                          neighbor =
                            connectivity.row_begin(additional),
                            end = connectivity.row_end(additional);
                          for (; neighbor!=end ; ++neighbor)
                            {
                              if (cell_partition[*neighbor] == part &&
                                  cell_partition_l2[*neighbor] ==
                                  size_info.n_active_cells)
                                {
                                  unsigned int this_index = 0;
                                  if (hp_bool == true)
                                    this_index = cell_active_fe_index.empty() ? 0 :
                                                 cell_active_fe_index[*neighbor];

                                  // Only add this cell if we need more macro
                                  // cells in the current block or if there is
                                  // a macro cell with the FE index that is
                                  // not yet fully populated
                                  if (missing_macros > 0 ||
                                      remaining_per_macro_cell[this_index] > 0)
                                    {
                                      cell_partition_l2[*neighbor] = partition_l2;
                                      neighbor_neighbor_list.push_back(*neighbor);
                                      if (hp_bool == true)
                                        renumbering_fe_index[this_index].
                                        push_back(*neighbor);
                                      partition_partition_list[counter] =
                                        *neighbor;
                                      counter++;
                                      partition_counter++;
                                      if (remaining_per_macro_cell[this_index]
                                          == 0 && missing_macros > 0)
                                        missing_macros--;
                                      remaining_per_macro_cell[this_index]++;
                                      if (remaining_per_macro_cell[this_index]
                                          == vectorization_length)
                                        {
                                          remaining_per_macro_cell[this_index] = 0;
                                        }
                                      if (missing_macros == 0)
                                        {
                                          filled = true;
                                          for (unsigned int fe_ind=0;
                                               fe_ind<max_fe_index; ++fe_ind)
                                            if (remaining_per_macro_cell[fe_ind]!=0)
                                              filled = false;
                                        }
                                      if (filled == true)
                                        break;
                                    }
                                }
                            }
                        }
                      if (hp_bool == true)
                        {
                          // set the renumbering according to their active FE
                          // index within one partition-partition which was
                          // implicitly assumed above
                          cell = counter - partition_counter;
                          for (unsigned int j=0; j<max_fe_index; j++)
                            {
                              for (unsigned int jj=0; jj<renumbering_fe_index[j].
                                   size(); jj++)
                                renumbering[cell++] =
                                  renumbering_fe_index[j][jj];
                              if (renumbering_fe_index[j].size()%vectorization_length != 0)
                                irregular_cells[renumbering_fe_index[j].size()/
                                                vectorization_length+
                                                n_macro_cells_before] =
                                                  renumbering_fe_index[j].size()%vectorization_length;
                              n_macro_cells_before += (renumbering_fe_index[j].
                                                       size()+vectorization_length-1)/
                                                      vectorization_length;
                              renumbering_fe_index[j].resize(0);
                            }
                        }
                      else
                        {
                          n_macro_cells_before += partition_counter/vectorization_length;
                          if (partition_counter%vectorization_length != 0)
                            {
                              irregular_cells[n_macro_cells_before] =
                                partition_counter%vectorization_length;
                              n_macro_cells_before++;
                            }
                        }
                    }
                    task_info.partition_color_blocks_data.
                    push_back(n_macro_cells_before);
                    partition_l2++;
                  }
                neighbor_list = neighbor_neighbor_list;
                neighbor_neighbor_list.resize(0);
              }
            task_info.partition_color_blocks_row_index[part+1] =
              task_info.partition_color_blocks_row_index[part] + partition_l2;
          }
      }

      if (size_info.boundary_cells_end>0)
        size_info.boundary_cells_end = task_info.partition_color_blocks_data
                                       [task_info.partition_color_blocks_row_index[1]];

      if (hp_bool == false)
        renumbering.swap(partition_partition_list);
      irregular_cells.resize(n_macro_cells_before);
      size_info.n_macro_cells = n_macro_cells_before;

      task_info.evens = (partition+1)/2;
      task_info.odds  = partition/2;
      task_info.n_blocked_workers =
        task_info.odds-(task_info.odds+task_info.evens+1)%2;
      task_info.n_workers = task_info.evens+task_info.odds-
                            task_info.n_blocked_workers;
      task_info.partition_evens.resize(partition);
      task_info.partition_odds.resize(partition);
      task_info.partition_n_blocked_workers.resize(partition);
      task_info.partition_n_workers.resize(partition);
      for (unsigned int part=0; part<partition; part++)
        {
          task_info.partition_evens[part] =
            (task_info.partition_color_blocks_row_index[part+1]-
             task_info.partition_color_blocks_row_index[part]+1)/2;
          task_info.partition_odds[part] =
            (task_info.partition_color_blocks_row_index[part+1]-
             task_info.partition_color_blocks_row_index[part])/2;
          task_info.partition_n_blocked_workers[part] =
            task_info.partition_odds[part]-(task_info.partition_odds[part]+
                                            task_info.partition_evens[part]+1)%2;
          task_info.partition_n_workers[part] =
            task_info.partition_evens[part]+task_info.partition_odds[part]-
            task_info.partition_n_blocked_workers[part];
        }
    }


    namespace internal
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
    }


    void
    DoFInfo::make_connectivity_graph
    (const SizeInfo                  &size_info,
     const TaskInfo                  &task_info,
     const std::vector<unsigned int> &renumbering,
     const std::vector<unsigned int> &irregular_cells,
     const bool                       do_blocking,
     CompressedSimpleSparsityPattern &connectivity) const
    {
      AssertDimension (row_starts.size()-1, size_info.n_active_cells);
      const unsigned int n_rows =
        (vector_partitioner->local_range().second-
         vector_partitioner->local_range().first)
        + vector_partitioner->ghost_indices().n_elements();
      const unsigned int n_blocks = (do_blocking == true) ?
                                    task_info.n_blocks : size_info.n_active_cells;

      // first determine row lengths
      std::vector<unsigned int> row_lengths(n_rows);
      unsigned int cell_start = 0, mcell_start = 0;
      std::vector<unsigned int> scratch;
      for (unsigned int block = 0; block < n_blocks; ++block)
        {
          // if we have the blocking variant (used in the coloring scheme), we
          // want to build a graph with the blocks with interaction with
          // remote MPI processes up front. in the non-blocking variant, we do
          // not do this here. TODO: unify this approach!!!
          if (do_blocking == true)
            {
              scratch.clear();
              for (unsigned int mcell=mcell_start; mcell<
                   std::min(mcell_start+task_info.block_size,
                            size_info.n_macro_cells);
                   ++mcell)
                {
                  unsigned int n_comp = (irregular_cells[mcell]>0)
                                        ?irregular_cells[mcell]:size_info.vectorization_length;
                  for (unsigned int cell = cell_start; cell < cell_start+n_comp;
                       ++cell)
                    scratch.insert(scratch.end(),
                                   begin_indices(renumbering[cell]),
                                   end_indices(renumbering[cell]));
                  cell_start += n_comp;
                }
              std::sort(scratch.begin(), scratch.end());
              const unsigned int n_unique =
                std::unique(scratch.begin(), scratch.end())-scratch.begin();
              for (unsigned int i=0; i<n_unique; ++i)
                row_lengths[scratch[i]]++;
              mcell_start += task_info.block_size;
            }
          else
            {
              scratch.clear();
              scratch.insert(scratch.end(),
                             begin_indices(block), end_indices(block));
              std::sort(scratch.begin(), scratch.end());
              const unsigned int n_unique =
                std::unique(scratch.begin(), scratch.end())-scratch.begin();
              for (unsigned int i=0; i<n_unique; ++i)
                row_lengths[scratch[i]]++;
            }
        }

      // disregard dofs that only sit on one cell
      for (unsigned int row=0; row<n_rows; ++row)
        if (row_lengths[row] == 1)
          row_lengths[row] = 0;

      SparsityPattern connectivity_dof (n_rows, n_blocks, row_lengths);
      cell_start = 0, mcell_start = 0;
      for (unsigned int block = 0; block < n_blocks; ++block)
        {
          // if we have the blocking variant (used in the coloring scheme), we
          // want to build a graph with the blocks with interaction with
          // remote MPI processes up front. in the non-blocking variant, we do
          // not do this here. TODO: unify this approach!!!
          if (do_blocking == true)
            {
              for (unsigned int mcell=mcell_start; mcell<
                   std::min(mcell_start+task_info.block_size,
                            size_info.n_macro_cells);
                   ++mcell)
                {
                  unsigned int n_comp = (irregular_cells[mcell]>0)
                                        ?irregular_cells[mcell]:size_info.vectorization_length;
                  for (unsigned int cell = cell_start; cell < cell_start+n_comp;
                       ++cell)
                    {
                      const unsigned int
                      *it = begin_indices (renumbering[cell]),
                       *end_cell = end_indices (renumbering[cell]);
                      for ( ; it != end_cell; ++it)
                        if (row_lengths[*it]>0)
                          connectivity_dof.add(*it, block);
                    }
                  cell_start += n_comp;
                }
              mcell_start += task_info.block_size;
            }
          else
            {
              const unsigned int
              *it = begin_indices (block),
               *end_cell = end_indices (block);
              for ( ; it != end_cell; ++it)
                if (row_lengths[*it]>0)
                  connectivity_dof.add(*it, block);
            }
        }
      connectivity_dof.compress();

      connectivity.reinit (n_blocks, n_blocks);
      internal::ordered_vector row_entries;
      cell_start = 0;
      mcell_start = 0;
      for (unsigned int block=0;  block < n_blocks; ++block)
        {
          row_entries.clear();

          if (do_blocking==true)
            {
              for (unsigned int mcell=mcell_start; mcell<
                   std::min(mcell_start+task_info.block_size,
                            size_info.n_macro_cells);
                   ++mcell)
                {
                  unsigned int n_comp = (irregular_cells[mcell]>0)
                                        ?irregular_cells[mcell]:size_info.vectorization_length;
                  for (unsigned int cell = cell_start; cell < cell_start+n_comp;
                       ++cell)
                    {
                      // apply renumbering when we do blocking
                      const unsigned int
                      *it = begin_indices (renumbering[cell]),
                       *end_cell = end_indices (renumbering[cell]);
                      for ( ; it != end_cell; ++it)
                        if (row_lengths[*it] > 0)
                          {
                            SparsityPattern::iterator sp = connectivity_dof.begin(*it);
                            // jump over diagonal for square patterns
                            if (connectivity_dof.n_rows()==connectivity_dof.n_cols())
                              ++sp;
                            row_entries.reserve (row_entries.size() + end_cell - it);
                            std::vector<types::global_dof_index>::iterator insert_pos = row_entries.begin();
                            for ( ; sp != connectivity_dof.end(*it); ++sp)
                              if (sp->column() >= block)
                                break;
                              else
                                row_entries.insert (sp->column(), insert_pos);
                          }
                    }
                  cell_start +=n_comp;
                }
              mcell_start += task_info.block_size;
            }
          else
            {
              const unsigned int *it = begin_indices (block),
                                  * end_cell = end_indices (block);
              for ( ; it != end_cell; ++it)
                if (row_lengths[*it] > 0)
                  {
                    SparsityPattern::iterator sp = connectivity_dof.begin(*it);
                    // jump over diagonal for square patterns
                    if (connectivity_dof.n_rows()==connectivity_dof.n_cols())
                      ++sp;
                    row_entries.reserve (row_entries.size() + end_cell - it);
                    std::vector<types::global_dof_index>::iterator insert_pos = row_entries.begin();
                    for ( ; sp != connectivity_dof.end(*it); ++sp)
                      if (sp->column() >= block)
                        break;
                      else
                        row_entries.insert (sp->column(), insert_pos);
                  }
            }
          connectivity.add_entries (block, row_entries.begin(), row_entries.end());
        }
      connectivity.symmetrize ();
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
      memory += (row_starts.capacity()*sizeof(std_cxx11::array<unsigned int,3>));
      memory += MemoryConsumption::memory_consumption (dof_indices);
      memory += MemoryConsumption::memory_consumption (row_starts_plain_indices);
      memory += MemoryConsumption::memory_consumption (plain_dof_indices);
      memory += MemoryConsumption::memory_consumption (constraint_indicator);
      memory += MemoryConsumption::memory_consumption (*vector_partitioner);
      return memory;
    }



    template <typename STREAM>
    void
    DoFInfo::print_memory_consumption (STREAM         &out,
                                       const SizeInfo &size_info) const
    {
      out << "       Memory row starts indices:    ";
      size_info.print_memory_statistics
      (out, (row_starts.capacity()*sizeof(std_cxx11::array<unsigned int, 3>)));
      out << "       Memory dof indices:           ";
      size_info.print_memory_statistics
      (out, MemoryConsumption::memory_consumption (dof_indices));
      out << "       Memory constraint indicators: ";
      size_info.print_memory_statistics
      (out, MemoryConsumption::memory_consumption (constraint_indicator));
      out << "       Memory plain indices:         ";
      size_info.print_memory_statistics
      (out, MemoryConsumption::memory_consumption (row_starts_plain_indices)+
       MemoryConsumption::memory_consumption (plain_dof_indices));
      out << "       Memory vector partitioner:    ";
      size_info.print_memory_statistics
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
