// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

#ifndef dealii_partitioner_templates_h
#define dealii_partitioner_templates_h

#include <deal.II/base/config.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <type_traits>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
#ifndef DOXYGEN

#  ifdef DEAL_II_WITH_MPI

    template <typename Number>
    void
    Partitioner::export_to_ghosted_array_start(
      const unsigned int             communication_channel,
      const ArrayView<const Number>& locally_owned_array,
      const ArrayView<Number>&       temporary_storage,
      const ArrayView<Number>&       ghost_array,
      std::vector<MPI_Request>&      requests) const
    {
      AssertDimension(temporary_storage.size(), n_import_indices());
      Assert(ghost_array.size() == n_ghost_indices()
               || ghost_array.size() == n_ghost_indices_in_larger_set,
             ExcGhostIndexArrayHasWrongSize(ghost_array.size(),
                                            n_ghost_indices(),
                                            n_ghost_indices_in_larger_set));

      const unsigned int n_import_targets = import_targets_data.size();
      const unsigned int n_ghost_targets  = ghost_targets_data.size();

      if(n_import_targets > 0)
        AssertDimension(locally_owned_array.size(), local_size());

      Assert(requests.size() == 0,
             ExcMessage("Another operation seems to still be running. "
                        "Call update_ghost_values_finish() first."));

      // Need to send and receive the data. Use non-blocking communication,
      // where it is usually less overhead to first initiate the receive and
      // then actually send the data
      requests.resize(n_import_targets + n_ghost_targets);

      // as a ghost array pointer, put the data at the end of the given ghost
      // array in case we want to fill only a subset of the ghosts so that we
      // can move data to the right position in a forward loop in the _finish
      // function.
      AssertIndexRange(n_ghost_indices(), n_ghost_indices_in_larger_set + 1);
      const bool use_larger_set
        = (n_ghost_indices_in_larger_set > n_ghost_indices()
           && ghost_array.size() == n_ghost_indices_in_larger_set);
      Number* ghost_array_ptr
        = use_larger_set ? ghost_array.data() + n_ghost_indices_in_larger_set
                             - n_ghost_indices() :
                           ghost_array.data();

      for(unsigned int i = 0; i < n_ghost_targets; i++)
        {
          // allow writing into ghost indices even though we are in a
          // const function
          const int ierr
            = MPI_Irecv(ghost_array_ptr,
                        ghost_targets_data[i].second * sizeof(Number),
                        MPI_BYTE,
                        ghost_targets_data[i].first,
                        ghost_targets_data[i].first + communication_channel,
                        communicator,
                        &requests[i]);
          AssertThrowMPI(ierr);
          ghost_array_ptr += ghost_targets()[i].second;
        }

      Number* temp_array_ptr = temporary_storage.data();
      for(unsigned int i = 0; i < n_import_targets; i++)
        {
          // copy the data to be sent to the import_data field
          std::vector<std::pair<unsigned int, unsigned int>>::const_iterator
            my_imports
            = import_indices_data.begin()
              + import_indices_chunks_by_rank_data[i],
            end_my_imports = import_indices_data.begin()
                             + import_indices_chunks_by_rank_data[i + 1];
          unsigned int index = 0;
          for(; my_imports != end_my_imports; ++my_imports)
            for(unsigned int j = my_imports->first; j < my_imports->second; j++)
              temp_array_ptr[index++] = locally_owned_array[j];
          AssertDimension(index, import_targets_data[i].second);

          // start the send operations
          const int ierr
            = MPI_Isend(temp_array_ptr,
                        import_targets_data[i].second * sizeof(Number),
                        MPI_BYTE,
                        import_targets_data[i].first,
                        my_pid + communication_channel,
                        communicator,
                        &requests[n_ghost_targets + i]);
          AssertThrowMPI(ierr);
          temp_array_ptr += import_targets_data[i].second;
        }
    }

    template <typename Number>
    void
    Partitioner::export_to_ghosted_array_finish(
      const ArrayView<Number>&  ghost_array,
      std::vector<MPI_Request>& requests) const
    {
      Assert(ghost_array.size() == n_ghost_indices()
               || ghost_array.size() == n_ghost_indices_in_larger_set,
             ExcGhostIndexArrayHasWrongSize(ghost_array.size(),
                                            n_ghost_indices(),
                                            n_ghost_indices_in_larger_set));

      // wait for both sends and receives to complete, even though only
      // receives are really necessary. this gives (much) better performance
      AssertDimension(ghost_targets().size() + import_targets().size(),
                      requests.size());
      if(requests.size() > 0)
        {
          const int ierr = MPI_Waitall(
            requests.size(), requests.data(), MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }
      requests.resize(0);

      // in case we only sent a subset of indices, we now need to move the data
      // to the correct positions and delete the old content
      if(n_ghost_indices_in_larger_set > n_ghost_indices()
         && ghost_array.size() == n_ghost_indices_in_larger_set)
        {
          unsigned int offset
            = n_ghost_indices_in_larger_set - n_ghost_indices();
          // must copy ghost data into extended ghost array
          for(std::vector<std::pair<unsigned int, unsigned int>>::const_iterator
                my_ghosts
              = ghost_indices_subset_data.begin();
              my_ghosts != ghost_indices_subset_data.end();
              ++my_ghosts)
            if(offset > my_ghosts->first)
              for(unsigned int j = my_ghosts->first; j < my_ghosts->second;
                  ++j, ++offset)
                {
                  ghost_array[j]      = ghost_array[offset];
                  ghost_array[offset] = Number();
                }
            else
              {
                AssertDimension(offset, my_ghosts->first);
                break;
              }
        }
    }

    template <typename Number>
    void
    Partitioner::import_from_ghosted_array_start(
      const VectorOperation::values vector_operation,
      const unsigned int            communication_channel,
      const ArrayView<Number>&      ghost_array,
      const ArrayView<Number>&      temporary_storage,
      std::vector<MPI_Request>&     requests) const
    {
      AssertDimension(temporary_storage.size(), n_import_indices());
      Assert(ghost_array.size() == n_ghost_indices()
               || ghost_array.size() == n_ghost_indices_in_larger_set,
             ExcGhostIndexArrayHasWrongSize(ghost_array.size(),
                                            n_ghost_indices(),
                                            n_ghost_indices_in_larger_set));

      (void) vector_operation;

      // nothing to do for insert (only need to zero ghost entries in
      // compress_finish()). in debug mode we want to check consistency of the
      // inserted data, therefore the communication is still initialized.
      // Having different code in debug and optimized mode is somewhat
      // dangerous, but it really saves communication so do it anyway
#    ifndef DEBUG
      if(vector_operation == VectorOperation::insert)
        return;
#    endif

      // nothing to do when we neither have import
      // nor ghost indices.
      if(n_ghost_indices() == 0 && n_import_indices() == 0)
        return;

      const unsigned int n_import_targets = import_targets_data.size();
      const unsigned int n_ghost_targets  = ghost_targets_data.size();

      Assert(requests.size() == 0,
             ExcMessage("Another compress operation seems to still be running. "
                        "Call compress_finish() first."));

      // Need to send and receive the data. Use non-blocking communication,
      // where it is generally less overhead to first initiate the receive and
      // then actually send the data

      // set channels in different range from update_ghost_values channels
      const unsigned int channel = communication_channel + 401;
      requests.resize(n_import_targets + n_ghost_targets);

      // initiate the receive operations
      Number* temp_array_ptr = temporary_storage.data();
      for(unsigned int i = 0; i < n_import_targets; i++)
        {
          AssertThrow(
            static_cast<std::size_t>(import_targets_data[i].second)
                * sizeof(Number)
              < static_cast<std::size_t>(std::numeric_limits<int>::max()),
            ExcMessage("Index overflow: Maximum message size in MPI is 2GB. "
                       "The number of ghost entries times the size of 'Number' "
                       "exceeds this value. This is not supported."));
          const int ierr
            = MPI_Irecv(temp_array_ptr,
                        import_targets_data[i].second * sizeof(Number),
                        MPI_BYTE,
                        import_targets_data[i].first,
                        import_targets_data[i].first + channel,
                        communicator,
                        &requests[i]);
          AssertThrowMPI(ierr);
          temp_array_ptr += import_targets_data[i].second;
        }

      // initiate the send operations

      // in case we want to import only from a subset of the ghosts we want to
      // move the data to send to the front of the array
      AssertIndexRange(n_ghost_indices(), n_ghost_indices_in_larger_set + 1);
      Number* ghost_array_ptr = ghost_array.data();
      for(unsigned int i = 0; i < n_ghost_targets; i++)
        {
          // in case we only sent a subset of indices, we now need to move the data
          // to the correct positions and delete the old content
          if(n_ghost_indices_in_larger_set > n_ghost_indices()
             && ghost_array.size() == n_ghost_indices_in_larger_set)
            {
              std::vector<std::pair<unsigned int, unsigned int>>::const_iterator
                my_ghosts
                = ghost_indices_subset_data.begin()
                  + ghost_indices_subset_chunks_by_rank_data[i],
                end_my_ghosts
                = ghost_indices_subset_data.begin()
                  + ghost_indices_subset_chunks_by_rank_data[i + 1];
              unsigned int offset = 0;
              for(; my_ghosts != end_my_ghosts; ++my_ghosts)
                if(ghost_array_ptr + offset
                   != ghost_array.data() + my_ghosts->first)
                  for(unsigned int j = my_ghosts->first; j < my_ghosts->second;
                      ++j, ++offset)
                    {
                      ghost_array_ptr[offset] = ghost_array[j];
                      ghost_array[j]          = Number();
                    }
                else
                  offset += my_ghosts->second - my_ghosts->first;
              AssertDimension(offset, ghost_targets_data[i].second);
            }

          AssertThrow(
            static_cast<std::size_t>(ghost_targets_data[i].second)
                * sizeof(Number)
              < static_cast<std::size_t>(std::numeric_limits<int>::max()),
            ExcMessage("Index overflow: Maximum message size in MPI is 2GB. "
                       "The number of ghost entries times the size of 'Number' "
                       "exceeds this value. This is not supported."));
          const int ierr
            = MPI_Isend(ghost_array_ptr,
                        ghost_targets_data[i].second * sizeof(Number),
                        MPI_BYTE,
                        ghost_targets_data[i].first,
                        this_mpi_process() + channel,
                        communicator,
                        &requests[n_import_targets + i]);
          AssertThrowMPI(ierr);

          ghost_array_ptr += ghost_targets_data[i].second;
        }
    }

    namespace internal
    {
      // In the import_from_ghosted_array_finish we need to invoke abs() also
      // on unsigned data types, which is ill-formed on newer C++
      // standards. To avoid this, we use std::abs on default types but
      // simply return the number on unsigned types
      template <typename Number>
      typename std::enable_if<
        !std::is_unsigned<Number>::value,
        typename numbers::NumberTraits<Number>::real_type>::type
      get_abs(const Number a)
      {
        return std::abs(a);
      }

      template <typename Number>
      typename std::enable_if<std::is_unsigned<Number>::value, Number>::type
      get_abs(const Number a)
      {
        return a;
      }

    } // namespace internal

    template <typename Number>
    void
    Partitioner::import_from_ghosted_array_finish(
      const VectorOperation::values  vector_operation,
      const ArrayView<const Number>& temporary_storage,
      const ArrayView<Number>&       locally_owned_array,
      const ArrayView<Number>&       ghost_array,
      std::vector<MPI_Request>&      requests) const
    {
      AssertDimension(temporary_storage.size(), n_import_indices());
      Assert(ghost_array.size() == n_ghost_indices()
               || ghost_array.size() == n_ghost_indices_in_larger_set,
             ExcGhostIndexArrayHasWrongSize(ghost_array.size(),
                                            n_ghost_indices(),
                                            n_ghost_indices_in_larger_set));

      // in optimized mode, no communication was started, so leave the
      // function directly (and only clear ghosts)
#    ifndef DEBUG
      if(vector_operation == VectorOperation::insert)
        {
          Assert(
            requests.empty(),
            ExcInternalError("Did not expect a non-empty communication "
                             "request when inserting. Check that the same "
                             "vector_operation argument was passed to "
                             "import_from_ghosted_array_start as is passed "
                             "to import_from_ghosted_array_finish."));
#      ifdef DEAL_II_WITH_CXX17
          if constexpr(std::is_trivial<Number>::value)
#      else
          if(std::is_trivial<Number>::value)
#      endif
            std::memset(
              ghost_array.data(), 0, sizeof(Number) * ghost_array.size());
          else
            std::fill(
              ghost_array.data(), ghost_array.data() + ghost_array.size(), 0);
          return;
        }
#    endif

      // nothing to do when we neither have import nor ghost indices.
      if(n_ghost_indices() == 0 && n_import_indices() == 0)
        return;

      const unsigned int n_import_targets = import_targets_data.size();
      const unsigned int n_ghost_targets  = ghost_targets_data.size();

      if(vector_operation != dealii::VectorOperation::insert)
        AssertDimension(n_ghost_targets + n_import_targets, requests.size());

      // first wait for the receive to complete
      if(requests.size() > 0 && n_import_targets > 0)
        {
          AssertDimension(locally_owned_array.size(), local_size());
          const int ierr = MPI_Waitall(
            n_import_targets, requests.data(), MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);

          const Number* read_position = temporary_storage.data();
          std::vector<std::pair<unsigned int, unsigned int>>::const_iterator
            my_imports
            = import_indices_data.begin();

          // If the operation is no insertion, add the imported data to the
          // local values. For insert, nothing is done here (but in debug mode
          // we assert that the specified value is either zero or matches with
          // the ones already present
          if(vector_operation != dealii::VectorOperation::insert)
            for(; my_imports != import_indices_data.end(); ++my_imports)
              for(unsigned int j = my_imports->first; j < my_imports->second;
                  j++)
                locally_owned_array[j] += *read_position++;
          else
            for(; my_imports != import_indices_data.end(); ++my_imports)
              for(unsigned int j = my_imports->first; j < my_imports->second;
                  j++, read_position++)
                // Below we use relatively large precision in units in the last place (ULP) as
                // this Assert can be easily triggered in p::d::SolutionTransfer.
                // The rationale is that during interpolation on two elements sharing
                // the face, values on this face obtained from each side might
                // be different due to additions being done in different order.
                Assert(*read_position == Number()
                         || internal::get_abs(locally_owned_array[j]
                                              - *read_position)
                              <= internal::get_abs(locally_owned_array[j]
                                                   + *read_position)
                                   * 100000.
                                   * std::numeric_limits<
                                       typename numbers::NumberTraits<
                                         Number>::real_type>::epsilon(),
                       typename LinearAlgebra::distributed::Vector<
                         Number>::ExcNonMatchingElements(*read_position,
                                                         locally_owned_array[j],
                                                         my_pid));
          AssertDimension(read_position - temporary_storage.data(),
                          n_import_indices());
        }

      // wait for the send operations to complete
      if(requests.size() > 0 && n_ghost_targets > 0)
        {
          const int ierr = MPI_Waitall(
            n_ghost_targets, &requests[n_import_targets], MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }
      else
        AssertDimension(n_ghost_indices(), 0);

      // clear the ghost array in case we did not yet do that in the _start
      // function
      if(ghost_array.size() > 0)
        {
          Assert(ghost_array.begin() != nullptr, ExcInternalError());
#    ifdef DEAL_II_WITH_CXX17
          if constexpr(std::is_trivial<Number>::value)
#    else
          if(std::is_trivial<Number>::value)
#    endif
            std::memset(
              ghost_array.data(), 0, sizeof(Number) * n_ghost_indices());
          else
            std::fill(
              ghost_array.data(), ghost_array.data() + n_ghost_indices(), 0);
        }

      // clear the compress requests
      requests.resize(0);
    }

#  endif // ifdef DEAL_II_WITH_MPI
#endif   // ifndef DOXYGEN

  } // end of namespace MPI

} // end of namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
