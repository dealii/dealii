// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/mpi_consensus_algorithms.h>
#include <deal.II/base/mpi_large_count.h>
#include <deal.II/base/mpi_tags.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <boost/serialization/utility.hpp>

// In this file, we use offsetof, which is a macro. When compiling
// with C++20 modules, this presents a problem because we wrap all of
// namespace std -- and then don't have access to macros. As a
// consequence, we really do need the following #include, even when
// building modules:
#include <cstddef> // Do not convert for module purposes
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <vector>

#if defined(DEAL_II_WITH_MPI)
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <mpi.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#endif


DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  IndexSet
  create_evenly_distributed_partitioning(
    const unsigned int            my_partition_id,
    const unsigned int            n_partitions,
    const types::global_dof_index total_size)
  {
    static_assert(std::is_same_v<types::global_dof_index, IndexSet::size_type>,
                  "IndexSet::size_type must match types::global_dof_index for "
                  "using this function");
    const unsigned int remain = total_size % n_partitions;

    const IndexSet::size_type min_size = total_size / n_partitions;

    const IndexSet::size_type begin =
      min_size * my_partition_id + std::min(my_partition_id, remain);
    const IndexSet::size_type end =
      min_size * (my_partition_id + 1) + std::min(my_partition_id + 1, remain);
    IndexSet result(total_size);
    result.add_range(begin, end);
    return result;
  }

  namespace MPI
  {
    MinMaxAvg
    min_max_avg(const double my_value, const MPI_Comm mpi_communicator)
    {
      MinMaxAvg result;
      min_max_avg(ArrayView<const double>(my_value),
                  ArrayView<MinMaxAvg>(result),
                  mpi_communicator);

      return result;
    }



    std::vector<MinMaxAvg>
    min_max_avg(const std::vector<double> &my_values,
                const MPI_Comm             mpi_communicator)
    {
      std::vector<MinMaxAvg> results(my_values.size());
      min_max_avg(my_values, results, mpi_communicator);

      return results;
    }



#ifdef DEAL_II_WITH_MPI
    unsigned int
    n_mpi_processes(const MPI_Comm mpi_communicator)
    {
      if (job_supports_mpi())
        {
          int       n_jobs = 1;
          const int ierr   = MPI_Comm_size(mpi_communicator, &n_jobs);
          AssertThrowMPI(ierr);
          return n_jobs;
        }
      else
        return 1;
    }


    unsigned int
    this_mpi_process(const MPI_Comm mpi_communicator)
    {
      if (job_supports_mpi())
        {
          int       rank = 0;
          const int ierr = MPI_Comm_rank(mpi_communicator, &rank);
          AssertThrowMPI(ierr);
          return rank;
        }
      else
        return 0;
    }



    std::vector<unsigned int>
    mpi_processes_within_communicator(const MPI_Comm comm_large,
                                      const MPI_Comm comm_small)
    {
      if (Utilities::MPI::job_supports_mpi() == false)
        return std::vector<unsigned int>{0};

      const unsigned int rank = Utilities::MPI::this_mpi_process(comm_large);
      const unsigned int size = Utilities::MPI::n_mpi_processes(comm_small);

      std::vector<unsigned int> ranks(size);
      const int                 ierr = MPI_Allgather(
        &rank, 1, MPI_UNSIGNED, ranks.data(), 1, MPI_UNSIGNED, comm_small);
      AssertThrowMPI(ierr);

      return ranks;
    }



    MPI_Comm
    duplicate_communicator(const MPI_Comm mpi_communicator)
    {
      MPI_Comm  new_communicator;
      const int ierr = MPI_Comm_dup(mpi_communicator, &new_communicator);
      AssertThrowMPI(ierr);
      return new_communicator;
    }



    void
    free_communicator(MPI_Comm mpi_communicator)
    {
      // MPI_Comm_free will set the argument to MPI_COMM_NULL automatically.
      const int ierr = MPI_Comm_free(&mpi_communicator);
      AssertThrowMPI(ierr);
    }



    std::vector<IndexSet>
    create_ascending_partitioning(
      const MPI_Comm                comm,
      const types::global_dof_index locally_owned_size)
    {
      static_assert(
        std::is_same_v<types::global_dof_index, IndexSet::size_type>,
        "IndexSet::size_type must match types::global_dof_index for "
        "using this function");
      const unsigned int                     n_proc = n_mpi_processes(comm);
      const std::vector<IndexSet::size_type> sizes =
        all_gather(comm, locally_owned_size);
      const auto total_size =
        std::accumulate(sizes.begin(), sizes.end(), IndexSet::size_type(0));

      std::vector<IndexSet> res(n_proc, IndexSet(total_size));

      IndexSet::size_type begin = 0;
      for (unsigned int i = 0; i < n_proc; ++i)
        {
          res[i].add_range(begin, begin + sizes[i]);
          begin = begin + sizes[i];
        }

      return res;
    }



    IndexSet
    create_evenly_distributed_partitioning(
      const MPI_Comm                comm,
      const types::global_dof_index total_size)
    {
      const unsigned int this_proc = this_mpi_process(comm);
      const unsigned int n_proc    = n_mpi_processes(comm);

      return Utilities::create_evenly_distributed_partitioning(this_proc,
                                                               n_proc,
                                                               total_size);
    }



    std::unique_ptr<MPI_Datatype, void (*)(MPI_Datatype *)>
    create_mpi_data_type_n_bytes(const std::size_t n_bytes)
    {
      MPI_Datatype result;
      int ierr = LargeCount::Type_contiguous_c(n_bytes, MPI_BYTE, &result);
      AssertThrowMPI(ierr);
      ierr = MPI_Type_commit(&result);
      AssertThrowMPI(ierr);

      if constexpr (running_in_debug_mode())
        {
          MPI_Count size64;
          ierr = MPI_Type_size_x(result, &size64);
          AssertThrowMPI(ierr);

          Assert(size64 == static_cast<MPI_Count>(n_bytes), ExcInternalError());
        }

      // Now put the new data type into a std::unique_ptr with a custom
      // deleter. We call the std::unique_ptr constructor that as first
      // argument takes a pointer (here, a pointer to a copy of the `result`
      // object, and as second argument a pointer-to-function, for which
      // we here use a lambda function without captures that acts as the
      // 'deleter' object: it calls `MPI_Type_free` and then deletes the
      // pointer. To avoid a compiler warning about a null this pointer
      // in the lambda (which don't make sense: the lambda doesn't store
      // anything), we create the deleter first.
      auto deleter = [](MPI_Datatype *p) {
        if (p != nullptr)
          {
            const int ierr = MPI_Type_free(p);
            AssertNothrow(ierr == MPI_SUCCESS, ExcMPI(ierr));
            delete p;
          }
      };

      return std::unique_ptr<MPI_Datatype, void (*)(MPI_Datatype *)>(
        new MPI_Datatype(result), deleter);
    }



    std::vector<unsigned int>
    compute_point_to_point_communication_pattern(
      const MPI_Comm                   mpi_comm,
      const std::vector<unsigned int> &destinations)
    {
      const unsigned int myid    = Utilities::MPI::this_mpi_process(mpi_comm);
      const unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_comm);

      if constexpr (running_in_debug_mode())
        {
          for (const unsigned int destination : destinations)
            AssertIndexRange(destination, n_procs);
        }

      // Have a little function that checks if destinations provided
      // to the current process are unique. The way it does this is
      // to create a sorted list of destinations and then walk through
      // the list and look at successive elements -- if we find the
      // same number twice, we know that the destinations were not
      // unique
      const bool my_destinations_are_unique = [destinations]() {
        if (destinations.empty())
          return true;
        else
          {
            std::vector<unsigned int> my_destinations = destinations;
            std::sort(my_destinations.begin(), my_destinations.end());
            return (std::adjacent_find(my_destinations.begin(),
                                       my_destinations.end()) ==
                    my_destinations.end());
          }
      }();

      // If all processes report that they have unique destinations,
      // then we can short-cut the process using a consensus algorithm (which
      // is implemented only for the case of unique destinations):
      if (Utilities::MPI::min((my_destinations_are_unique ? 1 : 0), mpi_comm) ==
          1)
        {
          return ConsensusAlgorithms::nbx<char, char>(
            destinations, {}, {}, {}, mpi_comm);
        }

      // So we need to run a different algorithm, specifically one that
      // requires more memory -- MPI_Reduce_scatter_block will require memory
      // proportional to the number of processes involved; that function is
      // available for MPI 2.2 or later:
      static CollectiveMutex      mutex;
      CollectiveMutex::ScopedLock lock(mutex, mpi_comm);

      const int mpi_tag =
        internal::Tags::compute_point_to_point_communication_pattern;

      // Calculate the number of messages to send to each process
      std::vector<unsigned int> dest_vector(n_procs);
      for (const auto &el : destinations)
        ++dest_vector[el];

      // Find how many processes will send to this one
      // by reducing with sum and then scattering the
      // results over all processes
      unsigned int n_recv_from;
      const int    ierr = MPI_Reduce_scatter_block(
        dest_vector.data(), &n_recv_from, 1, MPI_UNSIGNED, MPI_SUM, mpi_comm);

      AssertThrowMPI(ierr);

      // Send myid to every process in `destinations` vector...
      std::vector<MPI_Request> send_requests(destinations.size());
      for (const auto &el : destinations)
        {
          const int ierr =
            MPI_Isend(&myid,
                      1,
                      MPI_UNSIGNED,
                      el,
                      mpi_tag,
                      mpi_comm,
                      send_requests.data() + (&el - destinations.data()));
          AssertThrowMPI(ierr);
        }


      // Receive `n_recv_from` times from the processes
      // who communicate with this one. Store the obtained id's
      // in the resulting vector
      std::vector<unsigned int> origins(n_recv_from);
      for (auto &el : origins)
        {
          const int ierr = MPI_Recv(&el,
                                    1,
                                    MPI_UNSIGNED,
                                    MPI_ANY_SOURCE,
                                    mpi_tag,
                                    mpi_comm,
                                    MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);
        }

      if (destinations.size() > 0)
        {
          const int ierr = MPI_Waitall(destinations.size(),
                                       send_requests.data(),
                                       MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }

      return origins;
    }



    unsigned int
    compute_n_point_to_point_communications(
      const MPI_Comm                   mpi_comm,
      const std::vector<unsigned int> &destinations)
    {
      // Have a little function that checks if destinations provided
      // to the current process are unique:
      const bool my_destinations_are_unique = [destinations]() {
        std::vector<unsigned int> my_destinations = destinations;
        const unsigned int        n_destinations  = my_destinations.size();
        std::sort(my_destinations.begin(), my_destinations.end());
        my_destinations.erase(std::unique(my_destinations.begin(),
                                          my_destinations.end()),
                              my_destinations.end());
        return (my_destinations.size() == n_destinations);
      }();

      // If all processes report that they have unique destinations,
      // then we can short-cut the process using a consensus algorithm:

      if (Utilities::MPI::min((my_destinations_are_unique ? 1 : 0), mpi_comm) ==
          1)
        {
          return ConsensusAlgorithms::nbx<char, char>(
                   destinations, {}, {}, {}, mpi_comm)
            .size();
        }
      else
        {
          const unsigned int n_procs =
            Utilities::MPI::n_mpi_processes(mpi_comm);

          if constexpr (running_in_debug_mode())
            {
              for (const unsigned int destination : destinations)
                {
                  AssertIndexRange(destination, n_procs);
                  Assert(
                    destination != Utilities::MPI::this_mpi_process(mpi_comm),
                    ExcMessage(
                      "There is no point in communicating with ourselves."));
                }
            }

          // Calculate the number of messages to send to each process
          std::vector<unsigned int> dest_vector(n_procs);
          for (const auto &el : destinations)
            ++dest_vector[el];

          // Find out how many processes will send to this one
          // MPI_Reduce_scatter(_block) does exactly this
          unsigned int n_recv_from = 0;

          const int ierr = MPI_Reduce_scatter_block(dest_vector.data(),
                                                    &n_recv_from,
                                                    1,
                                                    MPI_UNSIGNED,
                                                    MPI_SUM,
                                                    mpi_comm);

          AssertThrowMPI(ierr);

          return n_recv_from;
        }
    }



    namespace
    {
      // custom MIP_Op for calculate_collective_mpi_min_max_avg
      void
      max_reduce(const void *in_lhs_,
                 void       *inout_rhs_,
                 int        *len,
                 MPI_Datatype *)
      {
        const MinMaxAvg *in_lhs    = static_cast<const MinMaxAvg *>(in_lhs_);
        MinMaxAvg       *inout_rhs = static_cast<MinMaxAvg *>(inout_rhs_);

        for (int i = 0; i < *len; ++i)
          {
            inout_rhs[i].sum += in_lhs[i].sum;
            if (inout_rhs[i].min > in_lhs[i].min)
              {
                inout_rhs[i].min       = in_lhs[i].min;
                inout_rhs[i].min_index = in_lhs[i].min_index;
              }
            else if (inout_rhs[i].min == in_lhs[i].min)
              {
                // choose lower cpu index when tied to make operator commutative
                if (inout_rhs[i].min_index > in_lhs[i].min_index)
                  inout_rhs[i].min_index = in_lhs[i].min_index;
              }

            if (inout_rhs[i].max < in_lhs[i].max)
              {
                inout_rhs[i].max       = in_lhs[i].max;
                inout_rhs[i].max_index = in_lhs[i].max_index;
              }
            else if (inout_rhs[i].max == in_lhs[i].max)
              {
                // choose lower cpu index when tied to make operator commutative
                if (inout_rhs[i].max_index > in_lhs[i].max_index)
                  inout_rhs[i].max_index = in_lhs[i].max_index;
              }
          }
      }
    } // namespace



    void
    min_max_avg(const ArrayView<const double> &my_values,
                const ArrayView<MinMaxAvg>    &result,
                const MPI_Comm                 mpi_communicator)
    {
      // If MPI was not started, we have a serial computation and cannot run
      // the other MPI commands
      if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 1)
        {
          for (unsigned int i = 0; i < my_values.size(); ++i)
            {
              result[i].sum       = my_values[i];
              result[i].avg       = my_values[i];
              result[i].min       = my_values[i];
              result[i].max       = my_values[i];
              result[i].min_index = 0;
              result[i].max_index = 0;
            }
          return;
        }

      /*
       * A custom MPI datatype handle describing the memory layout of the
       * MinMaxAvg struct. Initialized on first pass control reaches the
       * static variable. So hopefully not initialized too early.
       */
      static MPI_Datatype type = []() {
        MPI_Datatype type;

        int lengths[] = {3, 2, 1};

        MPI_Aint displacements[] = {0,
                                    offsetof(MinMaxAvg, min_index),
                                    offsetof(MinMaxAvg, avg)};

        MPI_Datatype types[] = {MPI_DOUBLE, MPI_INT, MPI_DOUBLE};

        int ierr =
          MPI_Type_create_struct(3, lengths, displacements, types, &type);
        AssertThrowMPI(ierr);

        ierr = MPI_Type_commit(&type);
        AssertThrowMPI(ierr);

        /* Ensure that we free the allocated datatype again at the end of
         * the program run just before we call MPI_Finalize():*/
        InitFinalize::signals.at_mpi_finalize.connect([type]() mutable {
          int ierr = MPI_Type_free(&type);
          AssertThrowMPI(ierr);
        });

        return type;
      }();

      /*
       * A custom MPI op handle for our max_reduce function.
       * Initialized on first pass control reaches the static variable. So
       * hopefully not initialized too early.
       */
      static MPI_Op op = []() {
        MPI_Op op;

        int ierr =
          MPI_Op_create(reinterpret_cast<MPI_User_function *>(&max_reduce),
                        static_cast<int>(true),
                        &op);
        AssertThrowMPI(ierr);

        /* Ensure that we free the allocated op again at the end of the
         * program run just before we call MPI_Finalize():*/
        InitFinalize::signals.at_mpi_finalize.connect([op]() mutable {
          int ierr = MPI_Op_free(&op);
          AssertThrowMPI(ierr);
        });

        return op;
      }();

      AssertDimension(Utilities::MPI::min(my_values.size(), mpi_communicator),
                      Utilities::MPI::max(my_values.size(), mpi_communicator));

      AssertDimension(my_values.size(), result.size());

      // To avoid uninitialized values on some MPI implementations, provide
      // result with a default value already...
      MinMaxAvg dummy = {0.,
                         std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::lowest(),
                         0,
                         0,
                         0.};

      for (auto &i : result)
        i = dummy;

      const unsigned int my_id =
        dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
      const unsigned int numproc =
        dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);

      std::vector<MinMaxAvg> in(my_values.size());

      for (unsigned int i = 0; i < my_values.size(); ++i)
        {
          in[i].sum = in[i].min = in[i].max = my_values[i];
          in[i].min_index = in[i].max_index = my_id;
        }

      int ierr = MPI_Allreduce(
        in.data(), result.data(), my_values.size(), type, op, mpi_communicator);
      AssertThrowMPI(ierr);

      for (auto &r : result)
        r.avg = r.sum / numproc;
    }


#else

    unsigned int
    n_mpi_processes(const MPI_Comm)
    {
      return 1;
    }



    unsigned int
    this_mpi_process(const MPI_Comm)
    {
      return 0;
    }



    std::vector<unsigned int>
    mpi_processes_within_communicator(const MPI_Comm, const MPI_Comm)
    {
      return std::vector<unsigned int>{0};
    }



    std::vector<IndexSet>
    create_ascending_partitioning(
      const MPI_Comm /*comm*/,
      const types::global_dof_index locally_owned_size)
    {
      return std::vector<IndexSet>(1, complete_index_set(locally_owned_size));
    }

    IndexSet
    create_evenly_distributed_partitioning(
      const MPI_Comm /*comm*/,
      const types::global_dof_index total_size)
    {
      return complete_index_set(total_size);
    }



    MPI_Comm
    duplicate_communicator(const MPI_Comm mpi_communicator)
    {
      return mpi_communicator;
    }



    void
    free_communicator(MPI_Comm /*mpi_communicator*/)
    {}



    void
    min_max_avg(const ArrayView<const double> &my_values,
                const ArrayView<MinMaxAvg>    &result,
                const MPI_Comm)
    {
      AssertDimension(my_values.size(), result.size());

      for (unsigned int i = 0; i < my_values.size(); ++i)
        {
          result[i].sum       = my_values[i];
          result[i].avg       = my_values[i];
          result[i].min       = my_values[i];
          result[i].max       = my_values[i];
          result[i].min_index = 0;
          result[i].max_index = 0;
        }
    }

#endif


    MPI_InitFinalize::MPI_InitFinalize(int               &argc,
                                       char            **&argv,
                                       const unsigned int max_num_threads)
      : InitFinalize(argc,
                     argv,
                     InitializeLibrary::MPI | InitializeLibrary::Kokkos |
                       InitializeLibrary::SLEPc | InitializeLibrary::PETSc |
                       InitializeLibrary::Zoltan | InitializeLibrary::P4EST,
                     max_num_threads)
    {}



    bool
    job_supports_mpi()
    {
#ifdef DEAL_II_WITH_MPI
      int       MPI_has_been_started = 0;
      const int ierr                 = MPI_Initialized(&MPI_has_been_started);
      AssertThrowMPI(ierr);

      return (MPI_has_been_started > 0);
#else
      return false;
#endif
    }



    namespace
    {
      /**
       * An internal namespace used for Utilities::MPI::compute_index_owner()
       * and for Utilities::MPI::Partitioner::set_ghost_indices().
       */
      namespace ComputeIndexOwner
      {
        class FlexibleIndexStorage
        {
        public:
          using index_type = unsigned int;
          static const index_type invalid_index_value =
            numbers::invalid_unsigned_int;

          FlexibleIndexStorage(const bool use_vector = true);

          void
          reinit(const bool        use_vector,
                 const bool        index_range_contiguous,
                 const std::size_t size);

          void
          fill(const std::size_t start,
               const std::size_t end,
               const index_type &value);

          index_type &
          operator[](const std::size_t index);

          index_type
          operator[](const std::size_t index) const;

          bool
          entry_has_been_set(const std::size_t index) const;

        private:
          bool                              use_vector;
          std::size_t                       size;
          std::vector<index_type>           data;
          std::map<std::size_t, index_type> data_map;
        };



        /**
         * Dictionary class with basic partitioning in terms of a single
         * interval of fixed size known to all MPI ranks for two-stage index
         * lookup.
         */
        struct Dictionary
        {
          /**
           * The minimum grain size for the intervals.
           *
           * We choose to limit the smallest size an interval for the
           * two-stage lookup can have with the following two conflicting
           * goals in mind: On the one hand, we do not want intervals in the
           * dictionary to become too short. For uneven distributions of
           * unknowns (some ranks with several thousands of unknowns, others
           * with none), the lookup DoFs -> dictionary then involves sending
           * from one MPI rank to many other MPI ranks holding dictionary
           * intervals, leading to an exceedingly high number of messages some
           * ranks have to send. Also, fewer longer intervals are generally
           * more efficient to look up. On the other hand, a range size too
           * large leads to opposite effect of many messages that come into a
           * particular dictionary owner in the lookup DoFs ->
           * dictionary. With the current setting, we get at most 64 messages
           * coming to a single MPI rank in case there is 1 dof per MPI rank,
           * which is reasonably low. At the same time, uneven distributions
           * up to factors of 4096 can be handled with at most 64 messages as
           * well.
           */
          static constexpr unsigned int range_minimum_grain_size = 64;

          /**
           * Factor that determines if an index set is sparse or not. An index
           * set if sparse if less than 25% of the indices are owned by any
           * process. If the index set is sparse, we switch the internal storage
           * from a fast storage (vector) to a memory-efficient storage (map).
           */
          static constexpr unsigned int sparsity_factor = 4;


          /**
           * Set up the dictionary by computing the partitioning from the
           * global size and sending the rank information on locally owned
           * ranges to the owner of the dictionary part.
           */
          Dictionary(const IndexSet &owned_indices, const MPI_Comm comm);

          /**
           * A vector with as many entries as there are dofs in the dictionary
           * of the current process, and each entry containing the rank of the
           * owner of that dof in the IndexSet `owned_indices`. This is
           * queried in the index lookup, so we keep an expanded list.
           */
          FlexibleIndexStorage actually_owning_ranks;

          /**
           * A sorted vector containing the MPI ranks appearing in
           * `actually_owning_ranks`.
           */
          std::vector<unsigned int> actually_owning_rank_list;

          /**
           * The number of unknowns in the dictionary for on each MPI rank
           * used for the index space splitting. For simplicity of index
           * lookup without additional communication, this number is the same
           * on all MPI ranks.
           */
          types::global_dof_index dofs_per_process;

          /**
           * The local range of the global index space that is represented in
           * the dictionary, computed from `dofs_per_process`, the current
           * MPI rank, and range_minimum_grain_size.
           */
          std::pair<types::global_dof_index, types::global_dof_index>
            local_range;

          /**
           * The actual size, computed as the minimum of dofs_per_process and
           * the possible end of the index space. Equivalent to
           * `local_range.second - local_range.first`.
           */
          types::global_dof_index locally_owned_size;

          /**
           * The global size of the index space.
           */
          types::global_dof_index size;

          /**
           * The number of ranks the `owned_indices` IndexSet is distributed
           * among.
           */
          unsigned int n_dict_procs_in_owned_indices;

          /**
           * A stride to distribute the work more evenly over MPI ranks in
           * case the grain size forces us to have fewer ranges than we have
           * processes.
           */
          unsigned int stride_small_size;

          /**
           * Translate a global dof index to the MPI rank in the dictionary
           * using `dofs_per_process`. We multiply by `stride_small_size` to
           * ensure a balance over the MPI ranks due to the grain size.
           */
          unsigned int
          dof_to_dict_rank(const types::global_dof_index i);

          /**
           * Given an MPI rank id of an arbitrary process, return the index
           * offset where the local range of that process begins.
           */
          types::global_dof_index
          get_index_offset(const unsigned int rank);

          /**
           * Given the rank in the owned indices from `actually_owning_ranks`,
           * this returns the index of the rank in the
           * `actually_owning_rank_list`.
           */
          unsigned int
          get_owning_rank_index(const unsigned int rank_in_owned_indices,
                                const unsigned int guess = 0);

        private:
          /**
           * Compute the partition from the global size of the index space and
           * the number of ranks.
           */
          void
          partition(const IndexSet &owned_indices, const MPI_Comm comm);
        };



        /**
         * Specialization of ConsensusAlgorithms::Process for the context of
         * Utilities::MPI::compute_index_owner() and
         * Utilities::MPI::Partitioner::set_ghost_indices() with additional
         * payload.
         */
        class ConsensusAlgorithmsPayload
        {
        public:
          using RequestType = std::vector<
            std::pair<types::global_dof_index, types::global_dof_index>>;
          using AnswerType = std::vector<unsigned int>;

          /**
           * Constructor.
           */
          ConsensusAlgorithmsPayload(const IndexSet &owned_indices,
                                     const IndexSet &indices_to_look_up,
                                     const MPI_Comm  comm,
                                     std::vector<unsigned int> &owning_ranks,
                                     const bool track_index_requesters = false);

          /**
           * The index space which describes the locally owned space.
           */
          const IndexSet &owned_indices;

          /**
           * The indices which are "ghosts" on a given rank and should be
           * looked up in terms of their owner rank from owned_indices.
           */
          const IndexSet &indices_to_look_up;

          /**
           * The underlying MPI communicator.
           */
          const MPI_Comm comm;

          /**
           * The present MPI rank.
           */
          const unsigned int my_rank;

          /**
           * The total number of ranks participating in the MPI communicator
           * `comm`.
           */
          const unsigned int n_procs;

          /**
           * Controls whether we should record a list of ranks who sent
           * requests to the present MPI process when looking up their remote
           * indices, and what those indices were. If true, it will be added
           * into `requesters` and can be queried by `get_requesters()`.
           */
          const bool track_index_requesters;

          /**
           * The result of the index owner computation: To each index
           * contained in `indices_to_look_up`, this vector contains the MPI
           * rank of the owner in `owned_indices`.
           */
          std::vector<unsigned int> &owning_ranks;

          /**
           * The dictionary handling the requests.
           */
          Dictionary dict;

          /**
           * Keeps track of the origin of the requests. The layout of the data
           * structure is as follows: The outermost vector has as many entries
           * as Dictionary::actually_owning_rank_list and represents the
           * information we should send back to the owners from the present
           * dictionary entry. The second vector then collects a list of MPI
           * ranks that have requested data, using the rank in the first pair
           * entry and a list of index ranges as the second entry.
           */
          std::vector<std::vector<
            std::pair<unsigned int,
                      std::vector<std::pair<types::global_dof_index,
                                            types::global_dof_index>>>>>
            requesters;

          /**
           * Array to collect the indices to look up (first vector) and their
           * local index among indices (second vector), sorted by the rank in
           * the dictionary.
           */
          std::map<unsigned int,
                   std::pair<std::vector<types::global_dof_index>,
                             std::vector<unsigned int>>>
            indices_to_look_up_by_dict_rank;

          /**
           * Return the recipients of requests.
           */
          std::vector<unsigned int>
          compute_targets();

          /**
           * The function that creates a request to another process.
           */
          std::vector<
            std::pair<types::global_dof_index, types::global_dof_index>>
          create_request(const unsigned int other_rank);

          /**
           * The function that answers a request from another process.
           */
          std::vector<unsigned int>
          answer_request(
            const unsigned int                                     other_rank,
            const std::vector<std::pair<types::global_dof_index,
                                        types::global_dof_index>> &buffer_recv);

          /**
           * The function that processes an answer from an MPI process we
           * have sent a request to.
           */
          void
          process_answer(const unsigned int               other_rank,
                         const std::vector<unsigned int> &recv_buffer);

          /**
           * Resolve the origin of the requests by sending the information
           * accumulated in terms of the dictionary owners during the run of
           * the consensus algorithm back to the owner in the original
           * IndexSet. This requires some point-to-point communication.
           *
           * @return Map of processes and associated sets of indices
           *         that are requested from the current rank. In
           *         other words, this function returns for each rank
           *         that has requested information about indices
           *         owned by the current which indices it has
           *         requested about; the values of the map are
           *         therefore all subsets of the owned set of
           *         indices.
           */
          std::map<unsigned int, IndexSet>
          get_requesters();

        private:
          /**
           * Stores the index request in the `requesters` field. Given the
           * rank of the owner, we start with a guess for the index at the
           * owner's site. This is because we typically might look up on the
           * same rank several times in a row, hence avoiding the binary
           * search in Dictionary::get_owning_rank_index()). Once we know the
           * index at the owner, we fill the vector entry with the rank of the
           * request. Here, we utilize the fact that requests are processed
           * rank-by-rank, so we can simply look at the end of the vector
           * whether there is already some data stored or not. Finally, we
           * build ranges, again using that the index list is sorted and we
           * therefore only need to append at the end.
           */
          void
          append_index_origin(
            const types::global_dof_index index_within_dictionary,
            const unsigned int            rank_of_request,
            const unsigned int            rank_of_owner,
            unsigned int                 &owner_index_guess);
        };

        /* ------------------------- inline functions ----------------------- */

        inline unsigned int
        Dictionary::dof_to_dict_rank(const types::global_dof_index i)
        {
          // note: this formula is also explicitly used in
          // get_index_offset(), so keep the two in sync
          return (i / dofs_per_process) * stride_small_size;
        }


        inline types::global_dof_index
        Dictionary::get_index_offset(const unsigned int rank)
        {
          return std::min(dofs_per_process *
                            static_cast<types::global_dof_index>(
                              (rank + stride_small_size - 1) /
                              stride_small_size),
                          size);
        }



        inline unsigned int
        Dictionary::get_owning_rank_index(
          const unsigned int rank_in_owned_indices,
          const unsigned int guess)
        {
          AssertIndexRange(guess, actually_owning_rank_list.size());
          if (actually_owning_rank_list[guess] == rank_in_owned_indices)
            return guess;
          else
            {
              auto it = std::lower_bound(actually_owning_rank_list.begin(),
                                         actually_owning_rank_list.end(),
                                         rank_in_owned_indices);
              Assert(it != actually_owning_rank_list.end(), ExcInternalError());
              Assert(*it == rank_in_owned_indices, ExcInternalError());
              return it - actually_owning_rank_list.begin();
            }
        }


        const FlexibleIndexStorage::index_type
          FlexibleIndexStorage::invalid_index_value;



        FlexibleIndexStorage::FlexibleIndexStorage(const bool use_vector)
          : use_vector(use_vector)
          , size(0)
        {}



        void
        FlexibleIndexStorage::reinit(const bool        use_vector,
                                     const bool        index_range_contiguous,
                                     const std::size_t size)
        {
          this->use_vector = use_vector;
          this->size       = size;

          data = {};
          data_map.clear();

          // in case we have contiguous indices, only fill the vector upon
          // first request in `fill`
          if (use_vector && !index_range_contiguous)
            data.resize(size, invalid_index_value);
        }



        void
        FlexibleIndexStorage::fill(
          const std::size_t                       start,
          const std::size_t                       end,
          const FlexibleIndexStorage::index_type &value)
        {
          AssertIndexRange(start, size);
          AssertIndexRange(end, size + 1);

          if (use_vector)
            {
              if (data.empty() && end > start)
                {
                  // in debug mode, we want to track whether we set all
                  // indices, so we first fill an invalid index and only later
                  // the actual ones, whereas we simply assign the given rank
                  // to the complete vector the first time we pass around in
                  // this function in release mode to avoid touching data
                  // unnecessarily (and overwrite the smaller pieces), as the
                  // locally owned part comes first
                  if constexpr (running_in_debug_mode())
                    {
                      data.resize(size, invalid_index_value);
                      std::fill(data.begin() + start,
                                data.begin() + end,
                                value);
                    }
                  else
                    {
                      data.resize(size, value);
                    }
                }
              else
                {
                  AssertDimension(data.size(), size);
                  std::fill(data.begin() + start, data.begin() + end, value);
                }
            }
          else
            {
              for (auto i = start; i < end; ++i)
                data_map[i] = value;
            }
        }



        FlexibleIndexStorage::index_type &
        FlexibleIndexStorage::operator[](const std::size_t index)
        {
          AssertIndexRange(index, size);

          if (use_vector)
            {
              AssertDimension(data.size(), size);
              return data[index];
            }
          else
            {
              return data_map.try_emplace(index, invalid_index_value)
                .first->second;
            }
        }



        inline bool
        FlexibleIndexStorage::entry_has_been_set(const std::size_t index) const
        {
          AssertIndexRange(index, size);

          if (use_vector)
            {
              if (data.empty())
                return false;

              AssertDimension(data.size(), size);
              return data[index] != invalid_index_value;
            }
          else
            return data_map.find(index) != data_map.end();
        }



        Dictionary::Dictionary(const IndexSet &owned_indices,
                               const MPI_Comm  comm)
        {
          // 1) set up the partition
          this->partition(owned_indices, comm);

          unsigned int my_rank = this_mpi_process(comm);

          types::global_dof_index dic_local_received = 0;
          std::map<unsigned int,
                   std::vector<std::pair<types::global_dof_index,
                                         types::global_dof_index>>>
            buffers;

          const auto owned_indices_size_actual =
            Utilities::MPI::sum(owned_indices.n_elements(), comm);

          actually_owning_ranks.reinit((owned_indices_size_actual *
                                        sparsity_factor) > owned_indices.size(),
                                       owned_indices_size_actual ==
                                         owned_indices.size(),
                                       locally_owned_size);

          // 2) collect relevant processes and process local dict entries
          for (auto interval = owned_indices.begin_intervals();
               interval != owned_indices.end_intervals();
               ++interval)
            {
              // Due to the granularity of the dictionary, the interval
              // might be split into several ranges of processor owner
              // ranks. Here, we process the interval by breaking into
              // smaller pieces in terms of the dictionary number.
              std::pair<types::global_dof_index, types::global_dof_index>
                index_range(*interval->begin(), interval->last() + 1);

              AssertThrow(index_range.second <= size, ExcInternalError());

              while (index_range.first != index_range.second)
                {
                  Assert(index_range.first < index_range.second,
                         ExcInternalError());

                  const unsigned int owner =
                    dof_to_dict_rank(index_range.first);

                  // this explicitly picks up the formula of
                  // dof_to_dict_rank, so the two places must be in sync
                  const types::global_dof_index next_index =
                    std::min(get_index_offset(owner + 1), index_range.second);

                  Assert(next_index > index_range.first, ExcInternalError());

                  if constexpr (running_in_debug_mode())
                    {
                      // make sure that the owner is the same on the current
                      // interval
                      for (types::global_dof_index i = index_range.first + 1;
                           i < next_index;
                           ++i)
                        AssertDimension(owner, dof_to_dict_rank(i));
                    }

                  // add the interval, either to the local range or into a
                  // buffer to be sent to another processor
                  if (owner == my_rank)
                    {
                      actually_owning_ranks.fill(index_range.first -
                                                   local_range.first,
                                                 next_index - local_range.first,
                                                 my_rank);
                      dic_local_received += next_index - index_range.first;
                      if (actually_owning_rank_list.empty())
                        actually_owning_rank_list.push_back(my_rank);
                    }
                  else
                    buffers[owner].emplace_back(index_range.first, next_index);

                  index_range.first = next_index;
                }
            }

#ifdef DEAL_II_WITH_MPI
          n_dict_procs_in_owned_indices = buffers.size();
          std::vector<MPI_Request> request;

          // Check if index set space is partitioned globally without gaps.
          if (owned_indices_size_actual == owned_indices.size())
            {
              // no gaps: setup is simple! Processes send their locally owned
              // indices to the dictionary. The dictionary stores the sending
              // rank for each index. The dictionary knows exactly
              // when it is set up when all indices it is responsible for
              // have been processed.

              request.reserve(n_dict_procs_in_owned_indices);

              // protect the following communication steps using a mutex:
              static CollectiveMutex      mutex;
              CollectiveMutex::ScopedLock lock(mutex, comm);

              const int mpi_tag =
                Utilities::MPI::internal::Tags::dictionary_reinit;


              // 3) send messages with local dofs to the right dict process
              for (const auto &rank_pair : buffers)
                {
                  request.push_back(MPI_Request());
                  const int ierr =
                    MPI_Isend(rank_pair.second.data(),
                              rank_pair.second.size() * 2,
                              Utilities::MPI::mpi_type_id_for_type<
                                types::global_dof_index>,
                              rank_pair.first,
                              mpi_tag,
                              comm,
                              &request.back());
                  AssertThrowMPI(ierr);
                }

              // 4) receive messages until all dofs in dict are processed
              while (this->locally_owned_size != dic_local_received)
                {
                  // wait for an incoming message
                  MPI_Status status;
                  int ierr = MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
                  AssertThrowMPI(ierr);

                  // retrieve size of incoming message
                  int number_amount;
                  ierr = MPI_Get_count(&status,
                                       Utilities::MPI::mpi_type_id_for_type<
                                         types::global_dof_index>,
                                       &number_amount);
                  AssertThrowMPI(ierr);

                  const auto other_rank = status.MPI_SOURCE;
                  actually_owning_rank_list.push_back(other_rank);

                  // receive message
                  Assert(number_amount % 2 == 0, ExcInternalError());
                  std::vector<
                    std::pair<types::global_dof_index, types::global_dof_index>>
                    buffer(number_amount / 2);
                  ierr = MPI_Recv(buffer.data(),
                                  number_amount,
                                  Utilities::MPI::mpi_type_id_for_type<
                                    types::global_dof_index>,
                                  status.MPI_SOURCE,
                                  status.MPI_TAG,
                                  comm,
                                  MPI_STATUS_IGNORE);
                  AssertThrowMPI(ierr);
                  // process message: loop over all intervals
                  for (auto interval : buffer)
                    {
                      if constexpr (library_build_mode ==
                                    LibraryBuildMode::debug)
                        {
                          for (types::global_dof_index i = interval.first;
                               i < interval.second;
                               i++)
                            Assert(actually_owning_ranks.entry_has_been_set(
                                     i - local_range.first) == false,
                                   ExcInternalError());
                          Assert(interval.first >= local_range.first &&
                                   interval.first < local_range.second,
                                 ExcInternalError());
                          Assert(interval.second > local_range.first &&
                                   interval.second <= local_range.second,
                                 ExcInternalError());
                        }

                      actually_owning_ranks.fill(interval.first -
                                                   local_range.first,
                                                 interval.second -
                                                   local_range.first,
                                                 other_rank);
                      dic_local_received += interval.second - interval.first;
                    }
                }
            }
          else
            {
              // with gap: use a ConsensusAlgorithm to determine when all
              // dictionaries have been set up.

              // 3/4) use a ConsensusAlgorithm to send messages with local
              // dofs to the right dict process

              using RequestType = std::vector<
                std::pair<types::global_dof_index, types::global_dof_index>>;

              ConsensusAlgorithms::selector<RequestType>(
                /* targets = */
                [&buffers]() {
                  std::vector<unsigned int> targets;
                  targets.reserve(buffers.size());
                  for (const auto &rank_pair : buffers)
                    targets.emplace_back(rank_pair.first);

                  return targets;
                }(),

                /* create_request = */
                [&buffers](const unsigned int target_rank) -> RequestType {
                  return buffers.at(target_rank);
                },

                /* process_request = */
                [&](const unsigned int source_rank,
                    const RequestType &request) -> void {
                  // process message: loop over all intervals
                  for (auto interval : request)
                    {
                      if constexpr (library_build_mode ==
                                    LibraryBuildMode::debug)
                        {
                          for (types::global_dof_index i = interval.first;
                               i < interval.second;
                               i++)
                            Assert(
                              actually_owning_ranks.entry_has_been_set(
                                i - local_range.first) == false,
                              ExcMessage(
                                "Multiple processes seem to own the same global index. "
                                "A possible reason is that the sets of locally owned "
                                "indices are not distinct."));
                          Assert(interval.first < interval.second,
                                 ExcInternalError());
                          Assert(
                            local_range.first <= interval.first &&
                              interval.second <= local_range.second,
                            ExcMessage(
                              "The specified interval is not handled by the current process."));
                        }
                      actually_owning_ranks.fill(interval.first -
                                                   local_range.first,
                                                 interval.second -
                                                   local_range.first,
                                                 source_rank);
                    }
                  actually_owning_rank_list.push_back(source_rank);
                },

                comm);
            }

          std::sort(actually_owning_rank_list.begin(),
                    actually_owning_rank_list.end());

          for (unsigned int i = 1; i < actually_owning_rank_list.size(); ++i)
            Assert(actually_owning_rank_list[i] >
                     actually_owning_rank_list[i - 1],
                   ExcInternalError());

          // 5) make sure that all messages have been sent
          if (request.size() > 0)
            {
              const int ierr = MPI_Waitall(request.size(),
                                           request.data(),
                                           MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
            }

#else
          Assert(buffers.empty(), ExcInternalError());
          (void)comm;
          (void)dic_local_received;
#endif
        }



        void
        Dictionary::partition(const IndexSet &owned_indices,
                              const MPI_Comm  comm)
        {
          const unsigned int n_procs = n_mpi_processes(comm);
          const unsigned int my_rank = this_mpi_process(comm);

          size = owned_indices.size();

          Assert(size > 0, ExcNotImplemented());

          dofs_per_process =
            std::max<types::global_dof_index>((size + n_procs - 1) / n_procs,
                                              range_minimum_grain_size);

          stride_small_size =
            std::max<unsigned int>(dofs_per_process * n_procs / size, 1);

          local_range.first  = get_index_offset(my_rank);
          local_range.second = get_index_offset(my_rank + 1);

          locally_owned_size = local_range.second - local_range.first;
        }


        ConsensusAlgorithmsPayload::ConsensusAlgorithmsPayload(
          const IndexSet            &owned_indices,
          const IndexSet            &indices_to_look_up,
          const MPI_Comm             comm,
          std::vector<unsigned int> &owning_ranks,
          const bool                 track_index_requesters)
          : owned_indices(owned_indices)
          , indices_to_look_up(indices_to_look_up)
          , comm(comm)
          , my_rank(this_mpi_process(comm))
          , n_procs(n_mpi_processes(comm))
          , track_index_requesters(track_index_requesters)
          , owning_ranks(owning_ranks)
          , dict(owned_indices, comm)
          , requesters(dict.actually_owning_rank_list.size())
        {}



        std::vector<unsigned int>
        ConsensusAlgorithmsPayload::compute_targets()
        {
          std::vector<unsigned int> targets;

          indices_to_look_up_by_dict_rank.clear();
          unsigned int index             = 0;
          unsigned int owner_index_guess = 0;
          for (auto i : indices_to_look_up)
            {
              unsigned int other_rank = dict.dof_to_dict_rank(i);
              if (other_rank == my_rank)
                {
                  owning_ranks[index] =
                    dict.actually_owning_ranks[i - dict.local_range.first];
                  if (track_index_requesters)
                    append_index_origin(i - dict.local_range.first,
                                        my_rank,
                                        owning_ranks[index],
                                        owner_index_guess);
                }
              else
                {
                  if (targets.empty() || targets.back() != other_rank)
                    targets.push_back(other_rank);
                  auto &indices = indices_to_look_up_by_dict_rank[other_rank];
                  indices.first.push_back(i);
                  indices.second.push_back(index);
                }
              ++index;
            }

          Assert(targets.size() == indices_to_look_up_by_dict_rank.size(),
                 ExcMessage("Size does not match!"));

          return targets;
        }



        std::vector<std::pair<types::global_dof_index, types::global_dof_index>>
        ConsensusAlgorithmsPayload::create_request(
          const unsigned int other_rank)
        {
          std::vector<
            std::pair<types::global_dof_index, types::global_dof_index>>
            send_buffer;

          // create index set and compress data to be sent
          auto &indices_i = indices_to_look_up_by_dict_rank[other_rank].first;
          IndexSet is(dict.size);
          is.add_indices(indices_i.begin(), indices_i.end());
          is.compress();

          for (auto interval = is.begin_intervals();
               interval != is.end_intervals();
               ++interval)
            send_buffer.emplace_back(*interval->begin(), interval->last() + 1);

          return send_buffer;
        }



        std::vector<unsigned int>
        ConsensusAlgorithmsPayload::answer_request(
          const unsigned int                                     other_rank,
          const std::vector<std::pair<types::global_dof_index,
                                      types::global_dof_index>> &buffer_recv)
        {
          std::vector<unsigned int> request_buffer;

          unsigned int owner_index_guess = 0;
          for (const auto &interval : buffer_recv)
            for (auto i = interval.first; i < interval.second; ++i)
              {
                const unsigned int actual_owner =
                  dict.actually_owning_ranks[i - dict.local_range.first];
                request_buffer.push_back(actual_owner);

                if (track_index_requesters)
                  append_index_origin(i - dict.local_range.first,
                                      other_rank,
                                      actual_owner,
                                      owner_index_guess);
              }

          return request_buffer;
        }



        void
        ConsensusAlgorithmsPayload::process_answer(
          const unsigned int               other_rank,
          const std::vector<unsigned int> &recv_buffer)
        {
          const auto &recv_indices =
            indices_to_look_up_by_dict_rank[other_rank].second;
          AssertDimension(recv_indices.size(), recv_buffer.size());
          for (unsigned int j = 0; j < recv_indices.size(); ++j)
            owning_ranks[recv_indices[j]] = recv_buffer[j];
        }



        std::map<unsigned int, IndexSet>
        ConsensusAlgorithmsPayload::get_requesters()
        {
          Assert(track_index_requesters,
                 ExcMessage("Must enable index range tracking in "
                            "constructor of ConsensusAlgorithmProcess"));

          std::map<unsigned int, dealii::IndexSet> requested_indices;

#ifdef DEAL_II_WITH_MPI

          static CollectiveMutex      mutex;
          CollectiveMutex::ScopedLock lock(mutex, comm);

          const int mpi_tag = Utilities::MPI::internal::Tags::
            consensus_algorithm_payload_get_requesters;

          // reserve enough slots for the requests ahead; depending on
          // whether the owning rank is one of the requesters or not, we
          // might have one less requests to execute, so fill the requests
          // on demand.
          std::vector<MPI_Request> send_requests;
          send_requests.reserve(requesters.size());

          // We use an integer vector for the data exchange. Since we send
          // data associated to intervals with different requesters, we will
          // need to send (a) the MPI rank of the requester, (b) the number
          // of intervals directed to this requester, and (c) a list of
          // intervals, i.e., two integers per interval. The number of items
          // sent in total can be deduced both via the MPI status message at
          // the receiver site as well as be counting the buckets from
          // different requesters.
          std::vector<std::vector<types::global_dof_index>> send_data(
            requesters.size());
          for (unsigned int i = 0; i < requesters.size(); ++i)
            {
              // special code for our own indices
              if (dict.actually_owning_rank_list[i] == my_rank)
                {
                  for (const auto &j : requesters[i])
                    {
                      const types::global_dof_index index_offset =
                        dict.get_index_offset(my_rank);
                      IndexSet &my_index_set = requested_indices[j.first];
                      my_index_set.set_size(owned_indices.size());
                      for (const auto &interval : j.second)
                        my_index_set.add_range(index_offset + interval.first,
                                               index_offset + interval.second);
                    }
                }
              else
                {
                  for (const auto &j : requesters[i])
                    {
                      send_data[i].push_back(j.first);
                      send_data[i].push_back(j.second.size());
                      for (const auto &interval : j.second)
                        {
                          send_data[i].push_back(interval.first);
                          send_data[i].push_back(interval.second);
                        }
                    }
                  send_requests.push_back(MPI_Request());
                  const int ierr =
                    MPI_Isend(send_data[i].data(),
                              send_data[i].size(),
                              Utilities::MPI::mpi_type_id_for_type<
                                types::global_dof_index>,
                              dict.actually_owning_rank_list[i],
                              mpi_tag,
                              comm,
                              &send_requests.back());
                  AssertThrowMPI(ierr);
                }
            }

          // receive the data
          for (unsigned int c = 0; c < dict.n_dict_procs_in_owned_indices; ++c)
            {
              // wait for an incoming message
              MPI_Status status;
              int ierr = MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
              AssertThrowMPI(ierr);

              // retrieve size of incoming message
              int number_amount;
              ierr = MPI_Get_count(
                &status,
                Utilities::MPI::mpi_type_id_for_type<types::global_dof_index>,
                &number_amount);
              AssertThrowMPI(ierr);

              // receive message
              Assert(number_amount % 2 == 0, ExcInternalError());
              std::vector<
                std::pair<types::global_dof_index, types::global_dof_index>>
                buffer(number_amount / 2);
              ierr = MPI_Recv(
                buffer.data(),
                number_amount,
                Utilities::MPI::mpi_type_id_for_type<types::global_dof_index>,
                status.MPI_SOURCE,
                status.MPI_TAG,
                comm,
                &status);
              AssertThrowMPI(ierr);

              // unpack the message and translate the dictionary-local
              // indices coming via MPI to the global index range
              const types::global_dof_index index_offset =
                dict.get_index_offset(status.MPI_SOURCE);
              unsigned int offset = 0;
              while (offset < buffer.size())
                {
                  AssertIndexRange(offset + buffer[offset].second,
                                   buffer.size());

                  IndexSet my_index_set(owned_indices.size());
                  for (unsigned int i = offset + 1;
                       i < offset + buffer[offset].second + 1;
                       ++i)
                    my_index_set.add_range(index_offset + buffer[i].first,
                                           index_offset + buffer[i].second);

                  // the underlying index set is able to merge ranges coming
                  // from different ranks due to the partitioning in the
                  // dictionary
                  IndexSet &index_set = requested_indices[buffer[offset].first];
                  if (index_set.size() == 0)
                    index_set.set_size(owned_indices.size());
                  index_set.add_indices(my_index_set);

                  offset += buffer[offset].second + 1;
                }
              AssertDimension(offset, buffer.size());
            }

          if (send_requests.size() > 0)
            {
              const auto ierr = MPI_Waitall(send_requests.size(),
                                            send_requests.data(),
                                            MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
            }


          if constexpr (running_in_debug_mode())
            {
              for (const auto &it : requested_indices)
                {
                  IndexSet copy_set = it.second;
                  copy_set.subtract_set(owned_indices);
                  Assert(copy_set.n_elements() == 0,
                         ExcInternalError(
                           "The indices requested from the current "
                           "MPI rank should be locally owned here!"));
                }
            }

#endif // DEAL_II_WITH_MPI

          return requested_indices;
        }



        void
        ConsensusAlgorithmsPayload::append_index_origin(
          const types::global_dof_index index_within_dict,
          const unsigned int            rank_of_request,
          const unsigned int            rank_of_owner,
          unsigned int                 &owner_index_guess)
        {
          // remember who requested which index. We want to use an
          // std::vector with simple addressing, via a good guess from the
          // preceding index, rather than std::map, because this is an inner
          // loop and it avoids the map lookup in every iteration
          owner_index_guess =
            dict.get_owning_rank_index(rank_of_owner, owner_index_guess);

          auto &request = requesters[owner_index_guess];
          if (request.empty() || request.back().first != rank_of_request)
            request.emplace_back(
              rank_of_request,
              std::vector<
                std::pair<types::global_dof_index, types::global_dof_index>>());

          auto &intervals = request.back().second;
          if (intervals.empty() || intervals.back().second != index_within_dict)
            intervals.emplace_back(index_within_dict, index_within_dict + 1);
          else
            ++intervals.back().second;
        }

      } // namespace ComputeIndexOwner
    }   // namespace



    std::vector<unsigned int>
    compute_index_owner(const IndexSet &owned_indices,
                        const IndexSet &indices_to_look_up,
                        const MPI_Comm  comm)
    {
      Assert(owned_indices.size() == indices_to_look_up.size(),
             ExcMessage("IndexSets have to have the same sizes."));

      Assert(
        owned_indices.size() == Utilities::MPI::max(owned_indices.size(), comm),
        ExcMessage("IndexSets have to have the same size on all processes."));

      std::vector<unsigned int> owning_ranks(indices_to_look_up.n_elements());

      // Step 1: setup dictionary
      // The input owned_indices can be partitioned arbitrarily. In the
      // dictionary, the index set is statically repartitioned among the
      // processes again and extended with information with the actual owner
      // of that the index.
      ComputeIndexOwner::ConsensusAlgorithmsPayload process(
        owned_indices,
        indices_to_look_up,
        comm,
        owning_ranks,
        /* keep track of requesters = */ false);

      // Step 2: read dictionary
      // Communicate with the process who owns the index in the static
      // partition (i.e. in the dictionary). This process returns the actual
      // owner of the index.
      using RequestType =
        ComputeIndexOwner::ConsensusAlgorithmsPayload::RequestType;
      using AnswerType =
        ComputeIndexOwner::ConsensusAlgorithmsPayload::AnswerType;
      ConsensusAlgorithms::selector<RequestType, AnswerType>(
        process.compute_targets(),
        [&process](const unsigned int other_rank) -> RequestType {
          return process.create_request(other_rank);
        },
        [&process](const unsigned int other_rank, const RequestType &r)
          -> AnswerType { return process.answer_request(other_rank, r); },
        [&process](const unsigned int other_rank, const AnswerType &a) -> void {
          process.process_answer(other_rank, a);
        },
        comm);

      return owning_ranks;
    }



    std::pair<std::vector<unsigned int>, std::map<unsigned int, IndexSet>>
    compute_index_owner_and_requesters(const IndexSet &owned_indices,
                                       const IndexSet &indices_to_look_up,
                                       const MPI_Comm &comm)
    {
      Assert(owned_indices.size() == indices_to_look_up.size(),
             ExcMessage("IndexSets have to have the same sizes."));

      Assert(
        owned_indices.size() == Utilities::MPI::max(owned_indices.size(), comm),
        ExcMessage("IndexSets have to have the same size on all processes."));

      std::vector<unsigned int> owning_ranks(indices_to_look_up.n_elements());

      // Step 1: setup dictionary
      // The input owned_indices can be partitioned arbitrarily. In the
      // dictionary, the index set is statically repartitioned among the
      // processes again and extended with information with the actual owner
      // of that the index.
      ComputeIndexOwner::ConsensusAlgorithmsPayload process(
        owned_indices, indices_to_look_up, comm, owning_ranks, true);

      // Step 2: read dictionary
      // Communicate with the process who owns the index in the static
      // partition (i.e. in the dictionary). This process returns the actual
      // owner of the index.
      using RequestType =
        ComputeIndexOwner::ConsensusAlgorithmsPayload::RequestType;
      using AnswerType =
        ComputeIndexOwner::ConsensusAlgorithmsPayload::AnswerType;
      ConsensusAlgorithms::selector<RequestType, AnswerType>(
        process.compute_targets(),
        [&process](const unsigned int other_rank) -> RequestType {
          return process.create_request(other_rank);
        },
        [&process](const unsigned int other_rank, const RequestType &r)
          -> AnswerType { return process.answer_request(other_rank, r); },
        [&process](const unsigned int other_rank, const AnswerType &a) -> void {
          process.process_answer(other_rank, a);
        },
        comm);

      return {owning_ranks, process.get_requesters()};
    }



    namespace internal
    {
      namespace CollectiveMutexImplementation
      {
        /**
         * Abort, should there be an exception being processed (see the error
         * message).
         */
        void
        check_exception()
        {
#ifdef DEAL_II_WITH_MPI
          if (std::uncaught_exceptions() > 0)
            {
              std::cerr
                << "---------------------------------------------------------\n"
                << "An exception was thrown inside a section of the program\n"
                << "guarded by a CollectiveMutex.\n"
                << "Because a CollectiveMutex guards critical communication\n"
                << "handling the exception would likely\n"
                << "deadlock because only the current process is aware of the\n"
                << "exception. To prevent this deadlock, the program will be\n"
                << "aborted.\n"
                << "---------------------------------------------------------"
                << std::endl;

              MPI_Abort(MPI_COMM_WORLD, 1);
            }
#endif
        }
      } // namespace CollectiveMutexImplementation
    }   // namespace internal



    CollectiveMutex::CollectiveMutex()
      : locked(false)
      , request(MPI_REQUEST_NULL)
    {
      InitFinalize::register_request(request);
    }



    CollectiveMutex::~CollectiveMutex()
    {
      // First check if this destructor is called during exception handling
      // if so, abort.
      internal::CollectiveMutexImplementation::check_exception();

      Assert(
        !locked,
        ExcMessage(
          "Error: MPI::CollectiveMutex is still locked while being destroyed!"));

      InitFinalize::unregister_request(request);
    }



    void
    CollectiveMutex::lock(const MPI_Comm comm)
    {
      Assert(
        !locked,
        ExcMessage(
          "Error: MPI::CollectiveMutex needs to be unlocked before lock()"));

#ifdef DEAL_II_WITH_MPI

      if (job_supports_mpi())
        {
          // TODO: For now, we implement this mutex with a blocking barrier in
          // the lock and unlock. It needs to be tested, if we can move to a
          // nonblocking barrier (code disabled below).

          const int ierr = MPI_Barrier(comm);
          AssertThrowMPI(ierr);

#  if 0
          // wait for non-blocking barrier to finish. This is a noop the
          // first time we lock().
          const int ierr = MPI_Wait(&request, MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);
#  else
          // nothing to do as blocking barrier already completed
#  endif
        }
#else
      (void)comm;
#endif

      locked = true;
    }



    void
    CollectiveMutex::unlock(const MPI_Comm comm)
    {
      // First check if this function is called during exception handling
      // if so, abort. This can happen if a ScopedLock is destroyed.
      internal::CollectiveMutexImplementation::check_exception();

      Assert(
        locked,
        ExcMessage(
          "Error: MPI::CollectiveMutex needs to be locked before unlock()"));

#ifdef DEAL_II_WITH_MPI

      if (job_supports_mpi())
        {
          // TODO: For now, we implement this mutex with a blocking barrier
          // in the lock and unlock. It needs to be tested, if we can move
          // to a nonblocking barrier (code disabled below):
#  if 0
      const int ierr = MPI_Ibarrier(comm, &request);
      AssertThrowMPI(ierr);
#  else
          const int ierr = MPI_Barrier(comm);
          AssertThrowMPI(ierr);
#  endif
        }
#else
      (void)comm;
#endif

      locked = false;
    }


#ifndef DOXYGEN
    // explicit instantiations

    // booleans aren't in MPI_SCALARS
    template bool
    reduce(const bool &,
           const MPI_Comm,
           const std::function<bool(const bool &, const bool &)> &,
           const unsigned int);

    template std::vector<bool>
    reduce(const std::vector<bool> &,
           const MPI_Comm,
           const std::function<std::vector<bool>(const std::vector<bool> &,
                                                 const std::vector<bool> &)> &,
           const unsigned int);

    template bool
    all_reduce(const bool &,
               const MPI_Comm,
               const std::function<bool(const bool &, const bool &)> &);

    template std::vector<bool>
    all_reduce(
      const std::vector<bool> &,
      const MPI_Comm,
      const std::function<std::vector<bool>(const std::vector<bool> &,
                                            const std::vector<bool> &)> &);

    // We need an explicit instantiation of this for the same reason as the
    // other types described in mpi.inst.in
    template void
    internal::all_reduce<bool>(const MPI_Op &,
                               const ArrayView<const bool> &,
                               const MPI_Comm,
                               const ArrayView<bool> &);


    template bool
    logical_or<bool>(const bool &, const MPI_Comm);


    template void
    logical_or<bool>(const ArrayView<const bool> &,
                     const MPI_Comm,
                     const ArrayView<bool> &);


    template std::vector<unsigned int>
    compute_set_union(const std::vector<unsigned int> &vec,
                      const MPI_Comm                   comm);


    template std::set<unsigned int>
    compute_set_union(const std::set<unsigned int> &set, const MPI_Comm comm);
#endif

#include "base/mpi.inst"
  } // end of namespace MPI
} // end of namespace Utilities

DEAL_II_NAMESPACE_CLOSE
