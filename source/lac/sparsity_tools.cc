// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2016 by the deal.II authors
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


#include <deal.II/base/exceptions.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <algorithm>
#include <functional>
#include <set>

#ifdef DEAL_II_WITH_MPI
#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#endif

#ifdef DEAL_II_WITH_METIS
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
extern "C"
{
#include <metis.h>
}
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#endif


DEAL_II_NAMESPACE_OPEN

namespace SparsityTools
{

  void partition (const SparsityPattern     &sparsity_pattern,
                  const unsigned int         n_partitions,
                  std::vector<unsigned int> &partition_indices)
  {
    Assert (sparsity_pattern.n_rows()==sparsity_pattern.n_cols(),
            ExcNotQuadratic());
    Assert (sparsity_pattern.is_compressed(),
            SparsityPattern::ExcNotCompressed());

    Assert (n_partitions > 0, ExcInvalidNumberOfPartitions(n_partitions));
    Assert (partition_indices.size() == sparsity_pattern.n_rows(),
            ExcInvalidArraySize (partition_indices.size(),
                                 sparsity_pattern.n_rows()));

    // check for an easy return
    if (n_partitions == 1 || (sparsity_pattern.n_rows()==1))
      {
        std::fill_n (partition_indices.begin(), partition_indices.size(), 0U);
        return;
      }

    // Make sure that METIS is actually
    // installed and detected
#ifndef DEAL_II_WITH_METIS
    (void)sparsity_pattern;
    AssertThrow (false, ExcMETISNotInstalled());
#else

    // generate the data structures for
    // METIS. Note that this is particularly
    // simple, since METIS wants exactly our
    // compressed row storage format. we only
    // have to set up a few auxiliary arrays
    idx_t
    n       = static_cast<signed int>(sparsity_pattern.n_rows()),
    ncon    = 1,                              // number of balancing constraints (should be >0)
    nparts  = static_cast<int>(n_partitions), // number of subdomains to create
    dummy;                                    // the numbers of edges cut by the
    // resulting partition

    // We can not partition n items into more than n parts. METIS will
    // generate non-sensical output (everything is owned by a single process)
    // and complain with a message (but won't return an error code!):
    // ***Cannot bisect a graph with 0 vertices!
    // ***You are trying to partition a graph into too many parts!
    nparts = std::min(n, nparts);

    // use default options for METIS
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions (options);

    // one more nuisance: we have to copy our own data to arrays that store
    // signed integers :-(
    std::vector<idx_t> int_rowstart(1);
    int_rowstart.reserve(sparsity_pattern.n_rows()+1);
    std::vector<idx_t> int_colnums;
    int_colnums.reserve(sparsity_pattern.n_nonzero_elements());
    for (SparsityPattern::size_type row=0; row<sparsity_pattern.n_rows(); ++row)
      {
        for (SparsityPattern::iterator col=sparsity_pattern.begin(row);
             col < sparsity_pattern.end(row); ++col)
          int_colnums.push_back(col->column());
        int_rowstart.push_back(int_colnums.size());
      }

    std::vector<idx_t> int_partition_indices (sparsity_pattern.n_rows());

    // Make use of METIS' error code.
    int ierr;

    // Select which type of partitioning to create

    // Use recursive if the number of partitions is less than or equal to 8
    if (nparts <= 8)
      ierr = METIS_PartGraphRecursive(&n, &ncon, &int_rowstart[0], &int_colnums[0],
                                      NULL, NULL, NULL,
                                      &nparts,NULL,NULL,&options[0],
                                      &dummy,&int_partition_indices[0]);

    // Otherwise use kway
    else
      ierr = METIS_PartGraphKway(&n, &ncon, &int_rowstart[0], &int_colnums[0],
                                 NULL, NULL, NULL,
                                 &nparts,NULL,NULL,&options[0],
                                 &dummy,&int_partition_indices[0]);

    // If metis returns normally, an error code METIS_OK=1 is returned from
    // the above functions (see metish.h)
    AssertThrow (ierr == 1, ExcMETISError (ierr));

    // now copy back generated indices into the output array
    std::copy (int_partition_indices.begin(),
               int_partition_indices.end(),
               partition_indices.begin());
#endif
  }


  namespace internal
  {
    /**
     * Given a connectivity graph and a list of indices (where
     * invalid_size_type indicates that a node has not been numbered yet),
     * pick a valid starting index among the as-yet unnumbered one.
     */
    DynamicSparsityPattern::size_type
    find_unnumbered_starting_index (const DynamicSparsityPattern &sparsity,
                                    const std::vector<DynamicSparsityPattern::size_type> &new_indices)
    {
      DynamicSparsityPattern::size_type starting_point   = numbers::invalid_size_type;
      DynamicSparsityPattern::size_type min_coordination = sparsity.n_rows();
      for (DynamicSparsityPattern::size_type row=0; row<sparsity.n_rows(); ++row)
        // look over all as-yet unnumbered indices
        if (new_indices[row] == numbers::invalid_size_type)
          {
            if (sparsity.row_length(row) < min_coordination)
              {
                min_coordination = sparsity.row_length(row);
                starting_point   = row;
              }
          }

      // now we still have to care for the case that no unnumbered dof has a
      // coordination number less than sparsity.n_rows(). this rather exotic
      // case only happens if we only have one cell, as far as I can see,
      // but there may be others as well.
      //
      // if that should be the case, we can chose an arbitrary dof as
      // starting point, e.g. the first unnumbered one
      if (starting_point == numbers::invalid_size_type)
        {
          for (DynamicSparsityPattern::size_type i=0; i<new_indices.size(); ++i)
            if (new_indices[i] == numbers::invalid_size_type)
              {
                starting_point = i;
                break;
              }

          Assert (starting_point != numbers::invalid_size_type,
                  ExcInternalError());
        }

      return starting_point;
    }
  }



  void
  reorder_Cuthill_McKee (const DynamicSparsityPattern                         &sparsity,
                         std::vector<DynamicSparsityPattern::size_type>       &new_indices,
                         const std::vector<DynamicSparsityPattern::size_type> &starting_indices)
  {
    Assert (sparsity.n_rows() == sparsity.n_cols(),
            ExcDimensionMismatch (sparsity.n_rows(), sparsity.n_cols()));
    Assert (sparsity.n_rows() == new_indices.size(),
            ExcDimensionMismatch (sparsity.n_rows(), new_indices.size()));
    Assert (starting_indices.size() <= sparsity.n_rows(),
            ExcMessage ("You can't specify more starting indices than there are rows"));
    Assert (sparsity.row_index_set().size() == 0 ||
            sparsity.row_index_set().size() == sparsity.n_rows(),
            ExcMessage("Only valid for sparsity patterns which store all rows."));
    for (SparsityPattern::size_type i=0; i<starting_indices.size(); ++i)
      Assert (starting_indices[i] < sparsity.n_rows(),
              ExcMessage ("Invalid starting index"));

    // store the indices of the dofs renumbered in the last round. Default to
    // starting points
    std::vector<DynamicSparsityPattern::size_type> last_round_dofs (starting_indices);

    // initialize the new_indices array with invalid values
    std::fill (new_indices.begin(), new_indices.end(),
               numbers::invalid_size_type);

    // delete disallowed elements
    for (DynamicSparsityPattern::size_type i=0; i<last_round_dofs.size(); ++i)
      if ((last_round_dofs[i]==numbers::invalid_size_type) ||
          (last_round_dofs[i]>=sparsity.n_rows()))
        last_round_dofs[i] = numbers::invalid_size_type;

    std::remove_if (last_round_dofs.begin(), last_round_dofs.end(),
                    std::bind2nd(std::equal_to<DynamicSparsityPattern::size_type>(),
                                 numbers::invalid_size_type));

    // now if no valid points remain: find dof with lowest coordination number
    if (last_round_dofs.empty())
      last_round_dofs
      .push_back (internal::find_unnumbered_starting_index (sparsity,
                                                            new_indices));

    // store next free dof index
    DynamicSparsityPattern::size_type next_free_number = 0;

    // enumerate the first round dofs
    for (DynamicSparsityPattern::size_type i=0; i!=last_round_dofs.size(); ++i)
      new_indices[last_round_dofs[i]] = next_free_number++;

    // now do as many steps as needed to renumber all dofs
    while (true)
      {
        // store the indices of the dofs to be renumbered in the next round
        std::vector<DynamicSparsityPattern::size_type> next_round_dofs;

        // find all neighbors of the dofs numbered in the last round
        for (DynamicSparsityPattern::size_type i=0; i<last_round_dofs.size(); ++i)
          for (DynamicSparsityPattern::iterator j=sparsity.begin(last_round_dofs[i]);
               j<sparsity.end(last_round_dofs[i]); ++j)
            next_round_dofs.push_back (j->column());

        // sort dof numbers
        std::sort (next_round_dofs.begin(), next_round_dofs.end());

        // delete multiple entries
        std::vector<DynamicSparsityPattern::size_type>::iterator end_sorted;
        end_sorted = std::unique (next_round_dofs.begin(), next_round_dofs.end());
        next_round_dofs.erase (end_sorted, next_round_dofs.end());

        // eliminate dofs which are already numbered
        for (int s=next_round_dofs.size()-1; s>=0; --s)
          if (new_indices[next_round_dofs[s]] != numbers::invalid_size_type)
            next_round_dofs.erase (next_round_dofs.begin() + s);

        // check whether there are any new dofs in the list. if there are
        // none, then we have completely numbered the current component of the
        // graph. check if there are as yet unnumbered components of the graph
        // that we would then have to do next
        if (next_round_dofs.empty())
          {
            if (std::find (new_indices.begin(), new_indices.end(),
                           numbers::invalid_size_type)
                ==
                new_indices.end())
              // no unnumbered indices, so we can leave now
              break;

            // otherwise find a valid starting point for the next component of
            // the graph and continue with numbering that one. we only do so
            // if no starting indices were provided by the user (see the
            // documentation of this function) so produce an error if we got
            // here and starting indices were given
            Assert (starting_indices.empty(),
                    ExcMessage ("The input graph appears to have more than one "
                                "component, but as stated in the documentation "
                                "we only want to reorder such graphs if no "
                                "starting indices are given. The function was "
                                "called with starting indices, however."))

            next_round_dofs
            .push_back (internal::find_unnumbered_starting_index (sparsity,
                                                                  new_indices));
          }



        // store for each coordination number the dofs with these coordination
        // number
        std::multimap<DynamicSparsityPattern::size_type, int> dofs_by_coordination;

        // find coordination number for each of these dofs
        for (std::vector<DynamicSparsityPattern::size_type>::iterator s=next_round_dofs.begin();
             s!=next_round_dofs.end(); ++s)
          {
            const DynamicSparsityPattern::size_type coordination = sparsity.row_length(*s);

            // insert this dof at its coordination number
            const std::pair<const DynamicSparsityPattern::size_type, int> new_entry (coordination, *s);
            dofs_by_coordination.insert (new_entry);
          }

        // assign new DoF numbers to the elements of the present front:
        std::multimap<DynamicSparsityPattern::size_type, int>::iterator i;
        for (i = dofs_by_coordination.begin(); i!=dofs_by_coordination.end(); ++i)
          new_indices[i->second] = next_free_number++;

        // after that: copy this round's dofs for the next round
        last_round_dofs = next_round_dofs;
      }

    // test for all indices numbered. this mostly tests whether the
    // front-marching-algorithm (which Cuthill-McKee actually is) has reached
    // all points.
    Assert ((std::find (new_indices.begin(), new_indices.end(), numbers::invalid_size_type)
             ==
             new_indices.end())
            &&
            (next_free_number == sparsity.n_rows()),
            ExcInternalError());
  }



  void
  reorder_Cuthill_McKee (const SparsityPattern                   &sparsity,
                         std::vector<SparsityPattern::size_type> &new_indices,
                         const std::vector<SparsityPattern::size_type> &starting_indices)
  {
    DynamicSparsityPattern dsp(sparsity.n_rows(), sparsity.n_cols());
    for (unsigned int row=0; row<sparsity.n_rows(); ++row)
      {
        for (SparsityPattern::iterator it=sparsity.begin(row); it!=sparsity.end(row)
             && it->is_valid_entry() ; ++it)
          dsp.add(row, it->column());
      }
    reorder_Cuthill_McKee(dsp, new_indices, starting_indices);
  }



  namespace internal
  {
    void
    reorder_hierarchical (const DynamicSparsityPattern                   &connectivity,
                          std::vector<DynamicSparsityPattern::size_type> &renumbering)
    {
      AssertDimension (connectivity.n_rows(), connectivity.n_cols());
      AssertDimension (connectivity.n_rows(), renumbering.size());
      Assert (connectivity.row_index_set().size() == 0 ||
              connectivity.row_index_set().size() == connectivity.n_rows(),
              ExcMessage("Only valid for sparsity patterns which store all rows."));

      std::vector<types::global_dof_index> touched_nodes(connectivity.n_rows(),
                                                         numbers::invalid_dof_index);
      std::vector<unsigned int> row_lengths(connectivity.n_rows());
      std::set<types::global_dof_index> current_neighbors;
      std::vector<std::vector<types::global_dof_index> > groups;

      // First collect the number of neighbors for each node. We use this
      // field to find next nodes with the minimum number of non-touched
      // neighbors in the field n_remaining_neighbors, so we will count down
      // on this field. We also cache the row lengths because we need this
      // data frequently and getting it from the sparsity pattern is more
      // expensive.
      for (types::global_dof_index row=0; row<connectivity.n_rows(); ++row)
        {
          row_lengths[row] = connectivity.row_length(row);
          Assert(row_lengths[row] > 0, ExcInternalError());
        }
      std::vector<unsigned int> n_remaining_neighbors(row_lengths);

      // This outer loop is typically traversed only once, unless the global
      // graph is not connected
      while (true)
        {
          // Find cell with the minimal number of neighbors (typically a
          // corner node when based on FEM meshes). If no cell is left, we are
          // done. Together with the outer while loop, this loop can possibly
          // be of quadratic complexity in the number of disconnected
          // partitions, i.e. up to connectivity.n_rows() in the worst case,
          // but that is not the usual use case of this loop and thus not
          // optimized for.
          std::pair<types::global_dof_index,types::global_dof_index> min_neighbors
          (numbers::invalid_dof_index, numbers::invalid_dof_index);
          for (types::global_dof_index i=0; i<touched_nodes.size(); ++i)
            if (touched_nodes[i] == numbers::invalid_dof_index)
              if (row_lengths[i] < min_neighbors.second)
                {
                  min_neighbors = std::make_pair(i, n_remaining_neighbors[i]);
                  if (n_remaining_neighbors[i] <= 1)
                    break;
                }
          if (min_neighbors.first == numbers::invalid_dof_index)
            break;

          Assert(min_neighbors.second > 0, ExcInternalError());

          current_neighbors.clear();
          current_neighbors.insert(min_neighbors.first);
          while (!current_neighbors.empty())
            {
              // Find node with minimum number of untouched neighbors among the
              // next set of possible neighbors
              min_neighbors = std::make_pair (numbers::invalid_dof_index,
                                              numbers::invalid_dof_index);
              for (std::set<types::global_dof_index>::iterator it=current_neighbors.begin();
                   it != current_neighbors.end(); ++it)
                {
                  Assert (touched_nodes[*it] == numbers::invalid_dof_index,
                          ExcInternalError());
                  if (n_remaining_neighbors[*it] < min_neighbors.second)
                    min_neighbors = std::make_pair(*it, n_remaining_neighbors[*it]);
                }

              // Among the set of nodes with the minimal number of neighbors,
              // choose the one with the largest number of touched neighbors,
              // i.e., the one with the largest row length
              const types::global_dof_index best_row_length = min_neighbors.second;
              for (std::set<types::global_dof_index>::iterator it=current_neighbors.begin();
                   it != current_neighbors.end(); ++it)
                if (n_remaining_neighbors[*it] == best_row_length)
                  if (row_lengths[*it] > min_neighbors.second)
                    min_neighbors = std::make_pair(*it, row_lengths[*it]);

              // Add the pivot and all direct neighbors of the pivot node not
              // yet touched to the list of new entries.
              groups.push_back(std::vector<types::global_dof_index>());
              std::vector<types::global_dof_index> &next_group = groups.back();

              next_group.push_back(min_neighbors.first);
              touched_nodes[min_neighbors.first] = groups.size()-1;
              for (DynamicSparsityPattern::iterator it
                   = connectivity.begin(min_neighbors.first);
                   it != connectivity.end(min_neighbors.first); ++it)
                if (touched_nodes[it->column()] == numbers::invalid_dof_index)
                  {
                    next_group.push_back(it->column());
                    touched_nodes[it->column()] = groups.size()-1;
                  }

              // Add all neighbors of the current list not yet touched to the
              // set of possible next pivots. The added node is no longer a
              // valid neighbor (here we assume symmetry of the
              // connectivity). Delete the entries of the current list from
              // the set of possible next pivots.
              for (unsigned int i=0; i<next_group.size(); ++i)
                {
                  for (DynamicSparsityPattern::iterator it
                       = connectivity.begin(next_group[i]);
                       it != connectivity.end(next_group[i]); ++it)
                    {
                      if (touched_nodes[it->column()] == numbers::invalid_dof_index)
                        current_neighbors.insert(it->column());
                      n_remaining_neighbors[it->column()]--;
                    }
                  current_neighbors.erase(next_group[i]);
                }
            }
        }

      // Sanity check: for all nodes, there should not be any neighbors left
      for (types::global_dof_index row=0; row<connectivity.n_rows(); ++row)
        Assert(n_remaining_neighbors[row] == 0, ExcInternalError());

      // If the number of groups is smaller than the number of nodes, we
      // continue by recursively calling this method
      if (groups.size() < connectivity.n_rows())
        {
          // Form the connectivity of the groups
          DynamicSparsityPattern connectivity_next(groups.size(),
                                                   groups.size());
          for (types::global_dof_index i=0; i<groups.size(); ++i)
            for (types::global_dof_index col=0; col<groups[i].size(); ++col)
              for (DynamicSparsityPattern::iterator it
                   = connectivity.begin(groups[i][col]);
                   it != connectivity.end(groups[i][col]); ++it)
                connectivity_next.add(i, touched_nodes[it->column()]);

          // Recursively call the reordering
          std::vector<types::global_dof_index> renumbering_next(groups.size());
          reorder_hierarchical(connectivity_next, renumbering_next);

          // Renumber the indices group by group according to the incoming
          // ordering for the groups
          for (types::global_dof_index i=0,count=0; i<groups.size(); ++i)
            for (types::global_dof_index col=0; col<groups[renumbering_next[i]].size(); ++col, ++count)
              renumbering[count] = groups[renumbering_next[i]][col];
        }
      else
        {
          // All groups should have size one and no more recursion is possible,
          // so use the numbering of the groups
          for (types::global_dof_index i=0,count=0; i<groups.size(); ++i)
            for (types::global_dof_index col=0; col<groups[i].size(); ++col, ++count)
              renumbering[count] = groups[i][col];
        }
    }
  }

  void
  reorder_hierarchical (const DynamicSparsityPattern                   &connectivity,
                        std::vector<DynamicSparsityPattern::size_type> &renumbering)
  {
    // the internal renumbering keeps the numbering the wrong way around (but
    // we cannot invert the numbering inside that method because it is used
    // recursively), so invert it here
    internal::reorder_hierarchical(connectivity, renumbering);
    renumbering = Utilities::invert_permutation(renumbering);
  }



#ifdef DEAL_II_WITH_MPI
  void distribute_sparsity_pattern
  (DynamicSparsityPattern                               &dsp,
   const std::vector<DynamicSparsityPattern::size_type> &rows_per_cpu,
   const MPI_Comm                                       &mpi_comm,
   const IndexSet                                       &myrange)
  {
    const unsigned int myid = Utilities::MPI::this_mpi_process(mpi_comm);
    std::vector<DynamicSparsityPattern::size_type> start_index(rows_per_cpu.size()+1);
    start_index[0]=0;
    for (DynamicSparsityPattern::size_type i=0; i<rows_per_cpu.size(); ++i)
      start_index[i+1]=start_index[i]+rows_per_cpu[i];

    typedef std::map<DynamicSparsityPattern::size_type,
            std::vector<DynamicSparsityPattern::size_type> >
            map_vec_t;

    map_vec_t send_data;

    {
      unsigned int dest_cpu=0;

      DynamicSparsityPattern::size_type n_local_rel_rows = myrange.n_elements();
      for (DynamicSparsityPattern::size_type row_idx=0; row_idx<n_local_rel_rows; ++row_idx)
        {
          DynamicSparsityPattern::size_type row=myrange.nth_index_in_set(row_idx);

          //calculate destination CPU
          while (row>=start_index[dest_cpu+1])
            ++dest_cpu;

          //skip myself
          if (dest_cpu==myid)
            {
              row_idx+=rows_per_cpu[myid]-1;
              continue;
            }

          DynamicSparsityPattern::size_type rlen = dsp.row_length(row);

          //skip empty lines
          if (!rlen)
            continue;

          //save entries
          std::vector<DynamicSparsityPattern::size_type> &dst = send_data[dest_cpu];

          dst.push_back(rlen); // number of entries
          dst.push_back(row); // row index
          for (DynamicSparsityPattern::size_type c=0; c<rlen; ++c)
            {
              //columns
              DynamicSparsityPattern::size_type column = dsp.column_number(row, c);
              dst.push_back(column);
            }
        }

    }

    unsigned int num_receive=0;
    {
      std::vector<unsigned int> send_to;
      send_to.reserve(send_data.size());
      for (map_vec_t::iterator it=send_data.begin(); it!=send_data.end(); ++it)
        send_to.push_back(it->first);

      num_receive =
        Utilities::MPI::
        compute_point_to_point_communication_pattern(mpi_comm, send_to).size();
    }

    std::vector<MPI_Request> requests(send_data.size());


    // send data
    {
      unsigned int idx=0;
      for (map_vec_t::iterator it=send_data.begin(); it!=send_data.end(); ++it, ++idx)
        MPI_Isend(&(it->second[0]),
                  it->second.size(),
                  DEAL_II_DOF_INDEX_MPI_TYPE,
                  it->first,
                  124,
                  mpi_comm,
                  &requests[idx]);
    }

    {
      //receive
      std::vector<DynamicSparsityPattern::size_type> recv_buf;
      for (unsigned int index=0; index<num_receive; ++index)
        {
          MPI_Status status;
          int len;
          MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpi_comm, &status);
          Assert (status.MPI_TAG==124, ExcInternalError());

          MPI_Get_count(&status, DEAL_II_DOF_INDEX_MPI_TYPE, &len);
          recv_buf.resize(len);
          MPI_Recv(&recv_buf[0], len, DEAL_II_DOF_INDEX_MPI_TYPE, status.MPI_SOURCE,
                   status.MPI_TAG, mpi_comm, &status);

          std::vector<DynamicSparsityPattern::size_type>::const_iterator ptr = recv_buf.begin();
          std::vector<DynamicSparsityPattern::size_type>::const_iterator end = recv_buf.end();
          while (ptr!=end)
            {
              DynamicSparsityPattern::size_type num=*(ptr++);
              Assert(ptr!=end, ExcInternalError());
              DynamicSparsityPattern::size_type row=*(ptr++);
              for (unsigned int c=0; c<num; ++c)
                {
                  Assert(ptr!=end, ExcInternalError());
                  dsp.add(row, *ptr);
                  ptr++;
                }
            }
          Assert(ptr==end, ExcInternalError());
        }
    }

    // complete all sends, so that we can safely destroy the buffers.
    if (requests.size())
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);

  }

  void distribute_sparsity_pattern(BlockDynamicSparsityPattern &dsp,
                                   const std::vector<IndexSet> &owned_set_per_cpu,
                                   const MPI_Comm              &mpi_comm,
                                   const IndexSet              &myrange)
  {
    const unsigned int myid = Utilities::MPI::this_mpi_process(mpi_comm);

    typedef std::map<BlockDynamicSparsityPattern::size_type,
            std::vector<BlockDynamicSparsityPattern::size_type> >
            map_vec_t;
    map_vec_t send_data;

    {
      unsigned int dest_cpu=0;

      BlockDynamicSparsityPattern::size_type n_local_rel_rows = myrange.n_elements();
      for (BlockDynamicSparsityPattern::size_type row_idx=0; row_idx<n_local_rel_rows; ++row_idx)
        {
          BlockDynamicSparsityPattern::size_type row=myrange.nth_index_in_set(row_idx);

          // calculate destination CPU, note that we start the search
          // at last destination cpu, because even if the owned ranges
          // are not contiguous, they hopefully consist of large blocks
          while (!owned_set_per_cpu[dest_cpu].is_element(row))
            {
              ++dest_cpu;
              if (dest_cpu==owned_set_per_cpu.size()) // wrap around
                dest_cpu=0;
            }

          //skip myself
          if (dest_cpu==myid)
            continue;

          BlockDynamicSparsityPattern::size_type rlen = dsp.row_length(row);

          //skip empty lines
          if (!rlen)
            continue;

          //save entries
          std::vector<BlockDynamicSparsityPattern::size_type> &dst = send_data[dest_cpu];

          dst.push_back(rlen); // number of entries
          dst.push_back(row); // row index
          for (BlockDynamicSparsityPattern::size_type c=0; c<rlen; ++c)
            {
              //columns
              BlockDynamicSparsityPattern::size_type column = dsp.column_number(row, c);
              dst.push_back(column);
            }
        }

    }

    unsigned int num_receive=0;
    {
      std::vector<unsigned int> send_to;
      send_to.reserve(send_data.size());
      for (map_vec_t::iterator it=send_data.begin(); it!=send_data.end(); ++it)
        send_to.push_back(it->first);

      num_receive =
        Utilities::MPI::
        compute_point_to_point_communication_pattern(mpi_comm, send_to).size();
    }

    std::vector<MPI_Request> requests(send_data.size());


    // send data
    {
      unsigned int idx=0;
      for (map_vec_t::iterator it=send_data.begin(); it!=send_data.end(); ++it, ++idx)
        MPI_Isend(&(it->second[0]),
                  it->second.size(),
                  DEAL_II_DOF_INDEX_MPI_TYPE,
                  it->first,
                  124,
                  mpi_comm,
                  &requests[idx]);
    }

    {
      //receive
      std::vector<BlockDynamicSparsityPattern::size_type> recv_buf;
      for (unsigned int index=0; index<num_receive; ++index)
        {
          MPI_Status status;
          int len;
          MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpi_comm, &status);
          Assert (status.MPI_TAG==124, ExcInternalError());

          MPI_Get_count(&status, DEAL_II_DOF_INDEX_MPI_TYPE, &len);
          recv_buf.resize(len);
          MPI_Recv(&recv_buf[0], len, DEAL_II_DOF_INDEX_MPI_TYPE, status.MPI_SOURCE,
                   status.MPI_TAG, mpi_comm, &status);

          std::vector<BlockDynamicSparsityPattern::size_type>::const_iterator ptr = recv_buf.begin();
          std::vector<BlockDynamicSparsityPattern::size_type>::const_iterator end = recv_buf.end();
          while (ptr!=end)
            {
              BlockDynamicSparsityPattern::size_type num=*(ptr++);
              Assert(ptr!=end, ExcInternalError());
              BlockDynamicSparsityPattern::size_type row=*(ptr++);
              for (unsigned int c=0; c<num; ++c)
                {
                  Assert(ptr!=end, ExcInternalError());
                  dsp.add(row, *ptr);
                  ptr++;
                }
            }
          Assert(ptr==end, ExcInternalError());
        }
    }

    // complete all sends, so that we can safely destroy the buffers.
    if (requests.size())
      MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
  }
#endif
}

DEAL_II_NAMESPACE_CLOSE
