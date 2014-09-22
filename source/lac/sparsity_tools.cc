// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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

#ifdef DEAL_II_WITH_MPI
#include <deal.II/base/utilities.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#endif

#ifdef DEAL_II_WITH_METIS
extern "C"
{
#include <metis.h>
}
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
    if (n_partitions == 1)
      {
        std::fill_n (partition_indices.begin(), partition_indices.size(), 0U);
        return;
      }

    // Make sure that METIS is actually
    // installed and detected
#ifndef DEAL_II_WITH_METIS
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

    // use default options for METIS
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions (options);

    // one more nuisance: we have to copy our
    // own data to arrays that store signed
    // integers :-(
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

    // Select which type of partitioning to
    // create

    // Use recursive if the number of
    // partitions is less than or equal to 8
    if (n_partitions <= 8)
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

    // If metis returns normally, an
    // error code METIS_OK=1 is
    // returned from the above
    // functions (see metish.h)
    AssertThrow (ierr == 1, ExcMETISError (ierr));

    // now copy back generated indices into the
    // output array
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
    SparsityPattern::size_type
    find_unnumbered_starting_index (const SparsityPattern     &sparsity,
                                    const std::vector<SparsityPattern::size_type> &new_indices)
    {
      {
        SparsityPattern::size_type starting_point   = numbers::invalid_size_type;
        SparsityPattern::size_type min_coordination = sparsity.n_rows();
        for (SparsityPattern::size_type row=0; row<sparsity.n_rows(); ++row)
          // look over all as-yet unnumbered indices
          if (new_indices[row] == numbers::invalid_size_type)
            {
              SparsityPattern::iterator j = sparsity.begin(row);

              // loop until we hit the end of this row's entries
              for ( ; j<sparsity.end(row); ++j)
                if (j->is_valid_entry() == false)
                  break;
              // post-condition after loop: coordination, i.e. the number of
              // entries in this row is now j-rowstart[row]
              if (static_cast<SparsityPattern::size_type>(j-sparsity.begin(row)) <
                  min_coordination)
                {
                  min_coordination = j-sparsity.begin(row);
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
            for (SparsityPattern::size_type i=0; i<new_indices.size(); ++i)
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
  }


  void
  reorder_Cuthill_McKee (const SparsityPattern                         &sparsity,
                         std::vector<SparsityPattern::size_type>       &new_indices,
                         const std::vector<SparsityPattern::size_type> &starting_indices)
  {
    Assert (sparsity.n_rows() == sparsity.n_cols(),
            ExcDimensionMismatch (sparsity.n_rows(), sparsity.n_cols()));
    Assert (sparsity.n_rows() == new_indices.size(),
            ExcDimensionMismatch (sparsity.n_rows(), new_indices.size()));
    Assert (starting_indices.size() <= sparsity.n_rows(),
            ExcMessage ("You can't specify more starting indices than there are rows"));
    for (SparsityPattern::size_type i=0; i<starting_indices.size(); ++i)
      Assert (starting_indices[i] < sparsity.n_rows(),
              ExcMessage ("Invalid starting index"));

    // store the indices of the dofs renumbered in the last round. Default to
    // starting points
    std::vector<SparsityPattern::size_type> last_round_dofs (starting_indices);

    // initialize the new_indices array with invalid values
    std::fill (new_indices.begin(), new_indices.end(),
               numbers::invalid_size_type);

    // delete disallowed elements
    for (SparsityPattern::size_type i=0; i<last_round_dofs.size(); ++i)
      if ((last_round_dofs[i]==numbers::invalid_size_type) ||
          (last_round_dofs[i]>=sparsity.n_rows()))
        last_round_dofs[i] = numbers::invalid_size_type;

    std::remove_if (last_round_dofs.begin(), last_round_dofs.end(),
                    std::bind2nd(std::equal_to<SparsityPattern::size_type>(),
                                 numbers::invalid_size_type));

    // now if no valid points remain: find dof with lowest coordination number
    if (last_round_dofs.empty())
      last_round_dofs
      .push_back (internal::find_unnumbered_starting_index (sparsity,
                                                            new_indices));

    // store next free dof index
    SparsityPattern::size_type next_free_number = 0;

    // enumerate the first round dofs
    for (SparsityPattern::size_type i=0; i!=last_round_dofs.size(); ++i)
      new_indices[last_round_dofs[i]] = next_free_number++;

    // now do as many steps as needed to
    // renumber all dofs
    while (true)
      {
        // store the indices of the dofs to be
        // renumbered in the next round
        std::vector<SparsityPattern::size_type> next_round_dofs;

        // find all neighbors of the
        // dofs numbered in the last
        // round
        for (SparsityPattern::size_type i=0; i<last_round_dofs.size(); ++i)
          for (SparsityPattern::iterator j=sparsity.begin(last_round_dofs[i]);
               j<sparsity.end(last_round_dofs[i]); ++j)
            if (j->is_valid_entry() == false)
              break;
            else
              next_round_dofs.push_back (j->column());

        // sort dof numbers
        std::sort (next_round_dofs.begin(), next_round_dofs.end());

        // delete multiple entries
        std::vector<SparsityPattern::size_type>::iterator end_sorted;
        end_sorted = std::unique (next_round_dofs.begin(), next_round_dofs.end());
        next_round_dofs.erase (end_sorted, next_round_dofs.end());

        // eliminate dofs which are
        // already numbered
        for (int s=next_round_dofs.size()-1; s>=0; --s)
          if (new_indices[next_round_dofs[s]] != numbers::invalid_size_type)
            next_round_dofs.erase (next_round_dofs.begin() + s);

        // check whether there are
        // any new dofs in the
        // list. if there are none,
        // then we have completely
        // numbered the current
        // component of the
        // graph. check if there are
        // as yet unnumbered
        // components of the graph
        // that we would then have to
        // do next
        if (next_round_dofs.empty())
          {
            if (std::find (new_indices.begin(), new_indices.end(),
                           numbers::invalid_size_type)
                ==
                new_indices.end())
              // no unnumbered
              // indices, so we can
              // leave now
              break;

            // otherwise find a valid
            // starting point for the
            // next component of the
            // graph and continue
            // with numbering that
            // one. we only do so if
            // no starting indices
            // were provided by the
            // user (see the
            // documentation of this
            // function) so produce
            // an error if we got
            // here and starting
            // indices were given
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



        // store for each coordination
        // number the dofs with these
        // coordination number
        std::multimap<SparsityPattern::size_type, int> dofs_by_coordination;

        // find coordination number for
        // each of these dofs
        for (std::vector<SparsityPattern::size_type>::iterator s=next_round_dofs.begin();
             s!=next_round_dofs.end(); ++s)
          {
            SparsityPattern::size_type coordination = 0;
            for (SparsityPattern::iterator j=sparsity.begin(*s);
                 j<sparsity.end(*s); ++j)
              if (j->is_valid_entry() == false)
                break;
              else
                ++coordination;

            // insert this dof at its
            // coordination number
            const std::pair<const SparsityPattern::size_type, int> new_entry (coordination, *s);
            dofs_by_coordination.insert (new_entry);
          }

        // assign new DoF numbers to
        // the elements of the present
        // front:
        std::multimap<SparsityPattern::size_type, int>::iterator i;
        for (i = dofs_by_coordination.begin(); i!=dofs_by_coordination.end(); ++i)
          new_indices[i->second] = next_free_number++;

        // after that: copy this round's
        // dofs for the next round
        last_round_dofs = next_round_dofs;
      }

    // test for all indices
    // numbered. this mostly tests
    // whether the
    // front-marching-algorithm (which
    // Cuthill-McKee actually is) has
    // reached all points.
    Assert ((std::find (new_indices.begin(), new_indices.end(), numbers::invalid_size_type)
             ==
             new_indices.end())
            &&
            (next_free_number == sparsity.n_rows()),
            ExcInternalError());
  }

#ifdef DEAL_II_WITH_MPI
  template <class CSP_t>
  void distribute_sparsity_pattern(CSP_t &csp,
                                   const std::vector<typename CSP_t::size_type> &rows_per_cpu,
                                   const MPI_Comm &mpi_comm,
                                   const IndexSet &myrange)
  {
    const unsigned int myid = Utilities::MPI::this_mpi_process(mpi_comm);
    std::vector<typename CSP_t::size_type> start_index(rows_per_cpu.size()+1);
    start_index[0]=0;
    for (typename CSP_t::size_type i=0; i<rows_per_cpu.size(); ++i)
      start_index[i+1]=start_index[i]+rows_per_cpu[i];

    typedef std::map<typename CSP_t::size_type, std::vector<typename CSP_t::size_type> > map_vec_t;

    map_vec_t send_data;

    {
      unsigned int dest_cpu=0;

      typename CSP_t::size_type n_local_rel_rows = myrange.n_elements();
      for (typename CSP_t::size_type row_idx=0; row_idx<n_local_rel_rows; ++row_idx)
        {
          typename CSP_t::size_type row=myrange.nth_index_in_set(row_idx);

          //calculate destination CPU
          while (row>=start_index[dest_cpu+1])
            ++dest_cpu;

          //skip myself
          if (dest_cpu==myid)
            {
              row_idx+=rows_per_cpu[myid]-1;
              continue;
            }

          typename CSP_t::size_type rlen = csp.row_length(row);

          //skip empty lines
          if (!rlen)
            continue;

          //save entries
          std::vector<typename CSP_t::size_type> &dst = send_data[dest_cpu];

          dst.push_back(rlen); // number of entries
          dst.push_back(row); // row index
          for (typename CSP_t::size_type c=0; c<rlen; ++c)
            {
              //columns
              typename CSP_t::size_type column = csp.column_number(row, c);
              dst.push_back(column);
            }
        }

    }

    unsigned int num_receive=0;
    {
      std::vector<unsigned int> send_to;
      send_to.reserve(send_data.size());
      for (typename map_vec_t::iterator it=send_data.begin(); it!=send_data.end(); ++it)
        send_to.push_back(it->first);

      num_receive =
        Utilities::MPI::
        compute_point_to_point_communication_pattern(mpi_comm, send_to).size();
    }

    std::vector<MPI_Request> requests(send_data.size());


    // send data
    {
      unsigned int idx=0;
      for (typename map_vec_t::iterator it=send_data.begin(); it!=send_data.end(); ++it, ++idx)
        MPI_Isend(&(it->second[0]),
                  it->second.size(),
                  DEAL_II_DOF_INDEX_MPI_TYPE,
                  it->first,
                  124,
                  mpi_comm,
                  &requests[idx]);
    }

//TODO: In the following, we read individual bytes and then reinterpret them
//    as typename CSP_t::size_type objects. this is error prone. use properly typed reads that
//    match the write above
    {
      //receive
      std::vector<typename CSP_t::size_type> recv_buf;
      for (unsigned int index=0; index<num_receive; ++index)
        {
          MPI_Status status;
          int len;
          MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpi_comm, &status);
          Assert (status.MPI_TAG==124, ExcInternalError());

          MPI_Get_count(&status, MPI_BYTE, &len);
          Assert( len%sizeof(typename CSP_t::size_type)==0, ExcInternalError());

          recv_buf.resize(len/sizeof(typename CSP_t::size_type));

          MPI_Recv(&recv_buf[0], len, MPI_BYTE, status.MPI_SOURCE,
                   status.MPI_TAG, mpi_comm, &status);

          typename std::vector<typename CSP_t::size_type>::const_iterator ptr = recv_buf.begin();
          typename std::vector<typename CSP_t::size_type>::const_iterator end = recv_buf.end();
          while (ptr+1<end)
            {
              typename CSP_t::size_type num=*(ptr++);
              typename CSP_t::size_type row=*(ptr++);
              for (unsigned int c=0; c<num; ++c)
                {
                  csp.add(row, *ptr);
                  ptr++;
                }
            }
          Assert(ptr==end, ExcInternalError());
        }
    }

    // complete all sends, so that we can safely destroy the buffers.
    MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);

  }

  template <class CSP_t>
  void distribute_sparsity_pattern(CSP_t &csp,
                                   const std::vector<IndexSet> &owned_set_per_cpu,
                                   const MPI_Comm &mpi_comm,
                                   const IndexSet &myrange)
  {
    const unsigned int myid = Utilities::MPI::this_mpi_process(mpi_comm);

    typedef std::map<typename CSP_t::size_type, std::vector<typename CSP_t::size_type> > map_vec_t;
    map_vec_t send_data;

    {
      unsigned int dest_cpu=0;

      typename CSP_t::size_type n_local_rel_rows = myrange.n_elements();
      for (typename CSP_t::size_type row_idx=0; row_idx<n_local_rel_rows; ++row_idx)
        {
          typename CSP_t::size_type row=myrange.nth_index_in_set(row_idx);

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

          typename CSP_t::size_type rlen = csp.row_length(row);

          //skip empty lines
          if (!rlen)
            continue;

          //save entries
          std::vector<typename CSP_t::size_type> &dst = send_data[dest_cpu];

          dst.push_back(rlen); // number of entries
          dst.push_back(row); // row index
          for (typename CSP_t::size_type c=0; c<rlen; ++c)
            {
              //columns
              typename CSP_t::size_type column = csp.column_number(row, c);
              dst.push_back(column);
            }
        }

    }

    unsigned int num_receive=0;
    {
      std::vector<unsigned int> send_to;
      send_to.reserve(send_data.size());
      for (typename map_vec_t::iterator it=send_data.begin(); it!=send_data.end(); ++it)
        send_to.push_back(it->first);

      num_receive =
        Utilities::MPI::
        compute_point_to_point_communication_pattern(mpi_comm, send_to).size();
    }

    std::vector<MPI_Request> requests(send_data.size());


    // send data
    {
      unsigned int idx=0;
      for (typename map_vec_t::iterator it=send_data.begin(); it!=send_data.end(); ++it, ++idx)
        MPI_Isend(&(it->second[0]),
                  it->second.size(),
                  DEAL_II_DOF_INDEX_MPI_TYPE,
                  it->first,
                  124,
                  mpi_comm,
                  &requests[idx]);
    }

//TODO: In the following, we read individual bytes and then reinterpret them
//    as typename CSP_t::size_type objects. this is error prone. use properly typed reads that
//    match the write above
    {
      //receive
      std::vector<typename CSP_t::size_type> recv_buf;
      for (unsigned int index=0; index<num_receive; ++index)
        {
          MPI_Status status;
          int len;
          MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpi_comm, &status);
          Assert (status.MPI_TAG==124, ExcInternalError());

          MPI_Get_count(&status, MPI_BYTE, &len);
          Assert( len%sizeof(typename CSP_t::size_type)==0, ExcInternalError());

          recv_buf.resize(len/sizeof(typename CSP_t::size_type));

          MPI_Recv(&recv_buf[0], len, MPI_BYTE, status.MPI_SOURCE,
                   status.MPI_TAG, mpi_comm, &status);

          typename std::vector<typename CSP_t::size_type>::const_iterator ptr = recv_buf.begin();
          typename std::vector<typename CSP_t::size_type>::const_iterator end = recv_buf.end();
          while (ptr+1<end)
            {
              typename CSP_t::size_type num=*(ptr++);
              typename CSP_t::size_type row=*(ptr++);
              for (unsigned int c=0; c<num; ++c)
                {
                  csp.add(row, *ptr);
                  ptr++;
                }
            }
          Assert(ptr==end, ExcInternalError());
        }
    }

    // complete all sends, so that we can safely destroy the buffers.
    MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
  }


#endif
}


//explicit instantiations

#define SPARSITY_FUNCTIONS(SparsityType) \
  template void SparsityTools::distribute_sparsity_pattern<SparsityType> (SparsityType & csp, \
      const std::vector<SparsityType::size_type> & rows_per_cpu,\
      const MPI_Comm & mpi_comm,\
      const IndexSet & myrange)

#ifdef DEAL_II_WITH_MPI
SPARSITY_FUNCTIONS(CompressedSparsityPattern);
SPARSITY_FUNCTIONS(CompressedSimpleSparsityPattern);

template void SparsityTools::distribute_sparsity_pattern
<BlockCompressedSimpleSparsityPattern>
(BlockCompressedSimpleSparsityPattern &csp,
 const std::vector<IndexSet> &owned_set_per_cpu,
 const MPI_Comm &mpi_comm,
 const IndexSet &myrange);

#endif

#undef SPARSITY_FUNCTIONS

DEAL_II_NAMESPACE_CLOSE
