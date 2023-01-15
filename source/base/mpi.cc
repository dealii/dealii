// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2022 by the deal.II authors
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


#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/mpi_large_count.h>
#include <deal.II/base/mpi_tags.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include <boost/serialization/utility.hpp>

#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <vector>

#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_WITH_MPI
#    include <deal.II/lac/trilinos_parallel_block_vector.h>
#    include <deal.II/lac/trilinos_vector.h>

#    include <Epetra_MpiComm.h>
#  endif
#endif

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/lac/petsc_block_vector.h>
#  include <deal.II/lac/petsc_vector.h>

#  include <petscsys.h>
#endif

#ifdef DEAL_II_WITH_SLEPC
#  include <deal.II/lac/slepc_solver.h>

#  include <slepcsys.h>
#endif

#ifdef DEAL_II_WITH_P4EST
#  include <p4est_bits.h>
#endif

#ifdef DEAL_II_TRILINOS_WITH_ZOLTAN
#  include <zoltan_cpp.h>
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
    static_assert(
      std::is_same<types::global_dof_index, IndexSet::size_type>::value,
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
#ifdef DEAL_II_WITH_MPI
    // Provide definitions of template variables for all valid instantiations.
    template const MPI_Datatype mpi_type_id_for_type<bool>;
    template const MPI_Datatype mpi_type_id_for_type<char>;
    template const MPI_Datatype mpi_type_id_for_type<signed char>;
    template const MPI_Datatype mpi_type_id_for_type<short>;
    template const MPI_Datatype mpi_type_id_for_type<int>;
    template const MPI_Datatype mpi_type_id_for_type<long int>;
    template const MPI_Datatype mpi_type_id_for_type<unsigned char>;
    template const MPI_Datatype mpi_type_id_for_type<unsigned short>;
    template const MPI_Datatype mpi_type_id_for_type<unsigned long int>;
    template const MPI_Datatype mpi_type_id_for_type<unsigned long long int>;
    template const MPI_Datatype mpi_type_id_for_type<float>;
    template const MPI_Datatype mpi_type_id_for_type<double>;
    template const MPI_Datatype mpi_type_id_for_type<long double>;
    template const MPI_Datatype mpi_type_id_for_type<std::complex<float>>;
    template const MPI_Datatype mpi_type_id_for_type<std::complex<double>>;
#endif


    MinMaxAvg
    min_max_avg(const double my_value, const MPI_Comm &mpi_communicator)
    {
      MinMaxAvg result;
      min_max_avg(ArrayView<const double>(my_value),
                  ArrayView<MinMaxAvg>(result),
                  mpi_communicator);

      return result;
    }



    std::vector<MinMaxAvg>
    min_max_avg(const std::vector<double> &my_values,
                const MPI_Comm &           mpi_communicator)
    {
      std::vector<MinMaxAvg> results(my_values.size());
      min_max_avg(my_values, results, mpi_communicator);

      return results;
    }



#ifdef DEAL_II_WITH_MPI
    unsigned int
    n_mpi_processes(const MPI_Comm &mpi_communicator)
    {
      int       n_jobs = 1;
      const int ierr   = MPI_Comm_size(mpi_communicator, &n_jobs);
      AssertThrowMPI(ierr);

      return n_jobs;
    }


    unsigned int
    this_mpi_process(const MPI_Comm &mpi_communicator)
    {
      int       rank = 0;
      const int ierr = MPI_Comm_rank(mpi_communicator, &rank);
      AssertThrowMPI(ierr);

      return rank;
    }



    const std::vector<unsigned int>
    mpi_processes_within_communicator(const MPI_Comm &comm_large,
                                      const MPI_Comm &comm_small)
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
    duplicate_communicator(const MPI_Comm &mpi_communicator)
    {
      MPI_Comm  new_communicator;
      const int ierr = MPI_Comm_dup(mpi_communicator, &new_communicator);
      AssertThrowMPI(ierr);
      return new_communicator;
    }



    void
    free_communicator(MPI_Comm &mpi_communicator)
    {
      // MPI_Comm_free will set the argument to MPI_COMM_NULL automatically.
      const int ierr = MPI_Comm_free(&mpi_communicator);
      AssertThrowMPI(ierr);
    }



    int
    create_group(const MPI_Comm & comm,
                 const MPI_Group &group,
                 const int        tag,
                 MPI_Comm *       new_comm)
    {
      const int ierr = MPI_Comm_create_group(comm, group, tag, new_comm);
      AssertThrowMPI(ierr);
      return ierr;
    }



    std::vector<IndexSet>
    create_ascending_partitioning(
      const MPI_Comm &              comm,
      const types::global_dof_index locally_owned_size)
    {
      static_assert(
        std::is_same<types::global_dof_index, IndexSet::size_type>::value,
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
      const MPI_Comm &              comm,
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

#  ifdef DEBUG
      MPI_Count size64;
      ierr = MPI_Type_size_x(result, &size64);
      AssertThrowMPI(ierr);

      Assert(size64 == static_cast<MPI_Count>(n_bytes), ExcInternalError());
#  endif

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
            (void)ierr;
            AssertNothrow(ierr == MPI_SUCCESS, ExcMPI(ierr));

            delete p;
          }
      };

      return std::unique_ptr<MPI_Datatype, void (*)(MPI_Datatype *)>(
        new MPI_Datatype(result), deleter);
    }



    std::vector<unsigned int>
    compute_point_to_point_communication_pattern(
      const MPI_Comm &                 mpi_comm,
      const std::vector<unsigned int> &destinations)
    {
      const unsigned int myid    = Utilities::MPI::this_mpi_process(mpi_comm);
      const unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_comm);
      (void)myid;
      (void)n_procs;

      for (const unsigned int destination : destinations)
        {
          (void)destination;
          AssertIndexRange(destination, n_procs);
        }


      // Have a little function that checks if destinations provided
      // to the current process are unique. The way it does this is
      // to create a sorted list of destinations and then walk through
      // the list and look at successive elements -- if we find the
      // same number twice, we know that the destinations were not
      // unique
      const bool my_destinations_are_unique = [destinations]() {
        if (destinations.size() == 0)
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
      const MPI_Comm &                 mpi_comm,
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

          for (const unsigned int destination : destinations)
            {
              (void)destination;
              AssertIndexRange(destination, n_procs);
              Assert(destination != Utilities::MPI::this_mpi_process(mpi_comm),
                     ExcMessage(
                       "There is no point in communicating with ourselves."));
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
                 void *      inout_rhs_,
                 int *       len,
                 MPI_Datatype *)
      {
        const MinMaxAvg *in_lhs    = static_cast<const MinMaxAvg *>(in_lhs_);
        MinMaxAvg *      inout_rhs = static_cast<MinMaxAvg *>(inout_rhs_);

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
                const ArrayView<MinMaxAvg> &   result,
                const MPI_Comm &               mpi_communicator)
    {
      // If MPI was not started, we have a serial computation and cannot run
      // the other MPI commands
      if (job_supports_mpi() == false ||
          Utilities::MPI::n_mpi_processes(mpi_communicator) <= 1)
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
        MPI_InitFinalize::signals.at_mpi_finalize.connect([type]() mutable {
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
        MPI_InitFinalize::signals.at_mpi_finalize.connect([op]() mutable {
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
    n_mpi_processes(const MPI_Comm &)
    {
      return 1;
    }



    unsigned int
    this_mpi_process(const MPI_Comm &)
    {
      return 0;
    }



    const std::vector<unsigned int>
    mpi_processes_within_communicator(const MPI_Comm &, const MPI_Comm &)
    {
      return std::vector<unsigned int>{0};
    }



    std::vector<IndexSet>
    create_ascending_partitioning(
      const MPI_Comm & /*comm*/,
      const types::global_dof_index locally_owned_size)
    {
      return std::vector<IndexSet>(1, complete_index_set(locally_owned_size));
    }

    IndexSet
    create_evenly_distributed_partitioning(
      const MPI_Comm & /*comm*/,
      const types::global_dof_index total_size)
    {
      return complete_index_set(total_size);
    }



    MPI_Comm
    duplicate_communicator(const MPI_Comm &mpi_communicator)
    {
      return mpi_communicator;
    }



    void
    free_communicator(MPI_Comm & /*mpi_communicator*/)
    {}



    void
    min_max_avg(const ArrayView<const double> &my_values,
                const ArrayView<MinMaxAvg> &   result,
                const MPI_Comm &)
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

    /* Force initialization of static struct: */
    MPI_InitFinalize::Signals MPI_InitFinalize::signals =
      MPI_InitFinalize::Signals();


    MPI_InitFinalize::MPI_InitFinalize(int &              argc,
                                       char **&           argv,
                                       const unsigned int max_num_threads)
    {
      static bool constructor_has_already_run = false;
      (void)constructor_has_already_run;
      Assert(constructor_has_already_run == false,
             ExcMessage("You can only create a single object of this class "
                        "in a program since it initializes the MPI system."));


      int ierr = 0;
#ifdef DEAL_II_WITH_MPI
      // if we have PETSc, we will initialize it and let it handle MPI.
      // Otherwise, we will do it.
      int MPI_has_been_started = 0;
      ierr                     = MPI_Initialized(&MPI_has_been_started);
      AssertThrowMPI(ierr);
      AssertThrow(MPI_has_been_started == 0,
                  ExcMessage("MPI error. You can only start MPI once!"));

      int provided;
      // this works like ierr = MPI_Init (&argc, &argv); but tells MPI that
      // we might use several threads but never call two MPI functions at the
      // same time. For an explanation see on why we do this see
      // http://www.open-mpi.org/community/lists/users/2010/03/12244.php
      int wanted = MPI_THREAD_SERIALIZED;
      ierr       = MPI_Init_thread(&argc, &argv, wanted, &provided);
      AssertThrowMPI(ierr);

      // disable for now because at least some implementations always return
      // MPI_THREAD_SINGLE.
      // Assert(max_num_threads==1 || provided != MPI_THREAD_SINGLE,
      //    ExcMessage("MPI reports that we are not allowed to use multiple
      //    threads."));
#else
      // make sure the compiler doesn't warn about these variables
      (void)argc;
      (void)argv;
      (void)ierr;
#endif

      // we are allowed to call MPI_Init ourselves and PETScInitialize will
      // detect this. This allows us to use MPI_Init_thread instead.
#ifdef DEAL_II_WITH_PETSC
#  ifdef DEAL_II_WITH_SLEPC
      // Initialize SLEPc (with PETSc):
      finalize_petscslepc = SlepcInitializeCalled ? false : true;
      ierr                = SlepcInitialize(&argc, &argv, nullptr, nullptr);
      AssertThrow(ierr == 0, SLEPcWrappers::SolverBase::ExcSLEPcError(ierr));
#  else
      // or just initialize PETSc alone:
      finalize_petscslepc = PetscInitializeCalled ? false : true;
      ierr                = PetscInitialize(&argc, &argv, nullptr, nullptr);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
#  endif

      // Disable PETSc exception handling. This just prints a large wall
      // of text that is not particularly helpful for what we do:
      PetscPopSignalHandler();
#endif

      // Initialize zoltan
#ifdef DEAL_II_TRILINOS_WITH_ZOLTAN
      float version;
      Zoltan_Initialize(argc, argv, &version);
#endif

#ifdef DEAL_II_WITH_P4EST
      // Initialize p4est and libsc components
#  if DEAL_II_P4EST_VERSION_GTE(2, 5, 0, 0)
      // This feature is broken in version 2.0.0 for calls to
      // MPI_Comm_create_group (see cburstedde/p4est#30).
      // Disabling it leads to more verbose p4est error messages
      // which should be fine.
      sc_init(MPI_COMM_WORLD, 0, 0, nullptr, SC_LP_SILENT);
#  endif
      p4est_init(nullptr, SC_LP_SILENT);
#endif

      constructor_has_already_run = true;


      // Now also see how many threads we'd like to run
      if (max_num_threads != numbers::invalid_unsigned_int)
        {
          // set maximum number of threads (also respecting the environment
          // variable that the called function evaluates) based on what the
          // user asked
          MultithreadInfo::set_thread_limit(max_num_threads);
        }
      else
        // user wants automatic choice
        {
#ifdef DEAL_II_WITH_MPI
          // we need to figure out how many MPI processes there are on the
          // current node, as well as how many CPU cores we have. for the
          // first task, check what get_hostname() returns and then do an
          // allgather so each processor gets the answer
          //
          // in calculating the length of the string, don't forget the
          // terminating \0 on C-style strings
          const std::string  hostname = Utilities::System::get_hostname();
          const unsigned int max_hostname_size =
            Utilities::MPI::max(hostname.size() + 1, MPI_COMM_WORLD);
          std::vector<char> hostname_array(max_hostname_size);
          std::copy(hostname.c_str(),
                    hostname.c_str() + hostname.size() + 1,
                    hostname_array.begin());

          std::vector<char> all_hostnames(max_hostname_size *
                                          MPI::n_mpi_processes(MPI_COMM_WORLD));
          const int         ierr = MPI_Allgather(hostname_array.data(),
                                         max_hostname_size,
                                         MPI_CHAR,
                                         all_hostnames.data(),
                                         max_hostname_size,
                                         MPI_CHAR,
                                         MPI_COMM_WORLD);
          AssertThrowMPI(ierr);

          // search how often our own hostname appears and the how-manyth
          // instance the current process represents
          unsigned int n_local_processes   = 0;
          unsigned int nth_process_on_host = 0;
          for (unsigned int i = 0; i < MPI::n_mpi_processes(MPI_COMM_WORLD);
               ++i)
            if (std::string(all_hostnames.data() + i * max_hostname_size) ==
                hostname)
              {
                ++n_local_processes;
                if (i <= MPI::this_mpi_process(MPI_COMM_WORLD))
                  ++nth_process_on_host;
              }
          Assert(nth_process_on_host > 0, ExcInternalError());


          // compute how many cores each process gets. if the number does not
          // divide evenly, then we get one more core if we are among the
          // first few processes
          //
          // if the number would be zero, round up to one since every process
          // needs to have at least one thread
          const unsigned int n_threads =
            std::max(MultithreadInfo::n_cores() / n_local_processes +
                       (nth_process_on_host <=
                            MultithreadInfo::n_cores() % n_local_processes ?
                          1 :
                          0),
                     1U);
#else
          const unsigned int n_threads = MultithreadInfo::n_cores();
#endif

          // finally set this number of threads
          MultithreadInfo::set_thread_limit(n_threads);
        }

      // As a final step call the at_mpi_init() signal handler.
      signals.at_mpi_init();
    }



    void
    MPI_InitFinalize::register_request(MPI_Request &request)
    {
      // insert if it is not in the set already:
      requests.insert(&request);
    }



    void
    MPI_InitFinalize::unregister_request(MPI_Request &request)
    {
      Assert(
        requests.find(&request) != requests.end(),
        ExcMessage(
          "You tried to call unregister_request() with an invalid request."));

      requests.erase(&request);
    }



    std::set<MPI_Request *> MPI_InitFinalize::requests;



    MPI_InitFinalize::~MPI_InitFinalize()
    {
      // First, call the at_mpi_finalize() signal handler.
      signals.at_mpi_finalize();

      // make memory pool release all PETSc/Trilinos/MPI-based vectors that
      // are no longer used at this point. this is relevant because the static
      // object destructors run for these vectors at the end of the program
      // would run after MPI_Finalize is called, leading to errors

#ifdef DEAL_II_WITH_MPI
      // Before exiting, wait for nonblocking communication to complete:
      for (auto request : requests)
        {
          const int ierr = MPI_Wait(request, MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);
        }

      // Start with deal.II MPI vectors and delete vectors from the pools:
      GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<double>>::release_unused_memory();
      GrowingVectorMemory<LinearAlgebra::distributed::BlockVector<double>>::
        release_unused_memory();
      GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<float>>::release_unused_memory();
      GrowingVectorMemory<LinearAlgebra::distributed::BlockVector<float>>::
        release_unused_memory();

      // Next with Trilinos:
#  ifdef DEAL_II_WITH_TRILINOS
      GrowingVectorMemory<
        TrilinosWrappers::MPI::Vector>::release_unused_memory();
      GrowingVectorMemory<
        TrilinosWrappers::MPI::BlockVector>::release_unused_memory();
#  endif
#endif


      // Now deal with PETSc (with or without MPI). Only delete the vectors if
      // finalize hasn't been called yet, otherwise this will lead to errors.
#ifdef DEAL_II_WITH_PETSC
      if (!PetscFinalizeCalled)
        {
          GrowingVectorMemory<
            PETScWrappers::MPI::Vector>::release_unused_memory();
          GrowingVectorMemory<
            PETScWrappers::MPI::BlockVector>::release_unused_memory();
        }
#  ifdef DEAL_II_WITH_SLEPC
      // and now end SLEPc with PETSc if we did so
      if (finalize_petscslepc)
        SlepcFinalize();
#  else
      // or just end PETSc if we did so
      if (finalize_petscslepc)
        PetscFinalize();
#  endif
#endif

#ifdef DEAL_II_WITH_P4EST
      // now end p4est and libsc
      // Note: p4est has no finalize function
      sc_finalize();
#endif


      // only MPI_Finalize if we are running with MPI. We also need to do this
      // when running PETSc, because we initialize MPI ourselves before
      // calling PetscInitialize
#ifdef DEAL_II_WITH_MPI
      if (job_supports_mpi() == true)
        {
#  if __cpp_lib_uncaught_exceptions >= 201411
          // std::uncaught_exception() is deprecated in c++17
          if (std::uncaught_exceptions() > 0)
#  else
          if (std::uncaught_exception() == true)
#  endif
            {
              // do not try to call MPI_Finalize to avoid a deadlock.
            }
          else
            {
              const int ierr = MPI_Finalize();
              (void)ierr;
              AssertNothrow(ierr == MPI_SUCCESS, dealii::ExcMPI(ierr));
            }
        }
#endif
    }



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



    std::vector<unsigned int>
    compute_index_owner(const IndexSet &owned_indices,
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
      internal::ComputeIndexOwner::ConsensusAlgorithmsPayload process(
        owned_indices, indices_to_look_up, comm, owning_ranks);

      // Step 2: read dictionary
      // Communicate with the process who owns the index in the static
      // partition (i.e. in the dictionary). This process returns the actual
      // owner of the index.
      ConsensusAlgorithms::Selector<
        std::vector<
          std::pair<types::global_dof_index, types::global_dof_index>>,
        std::vector<unsigned int>>
        consensus_algorithm;
      consensus_algorithm.run(process, comm);

      return owning_ranks;
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
#  if __cpp_lib_uncaught_exceptions >= 201411
          // std::uncaught_exception() is deprecated in c++17
          if (std::uncaught_exceptions() != 0)
#  else
          if (std::uncaught_exception() == true)
#  endif
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
      Utilities::MPI::MPI_InitFinalize::register_request(request);
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

      Utilities::MPI::MPI_InitFinalize::unregister_request(request);
    }



    void
    CollectiveMutex::lock(const MPI_Comm &comm)
    {
      (void)comm;

      Assert(
        !locked,
        ExcMessage(
          "Error: MPI::CollectiveMutex needs to be unlocked before lock()"));

#ifdef DEAL_II_WITH_MPI

      // TODO: For now, we implement this mutex with a blocking barrier
      // in the lock and unlock. It needs to be tested, if we can move
      // to a nonblocking barrier (code disabled below).

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
#endif

      locked = true;
    }



    void
    CollectiveMutex::unlock(const MPI_Comm &comm)
    {
      (void)comm;

      // First check if this function is called during exception handling
      // if so, abort. This can happen if a ScopedLock is destroyed.
      internal::CollectiveMutexImplementation::check_exception();

      Assert(
        locked,
        ExcMessage(
          "Error: MPI::CollectiveMutex needs to be locked before unlock()"));

#ifdef DEAL_II_WITH_MPI

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
#endif

      locked = false;
    }


#ifndef DOXYGEN
    // explicit instantiations

    // booleans aren't in MPI_SCALARS
    template bool
    reduce(const bool &,
           const MPI_Comm &,
           const std::function<bool(const bool &, const bool &)> &,
           const unsigned int);

    template std::vector<bool>
    reduce(const std::vector<bool> &,
           const MPI_Comm &,
           const std::function<std::vector<bool>(const std::vector<bool> &,
                                                 const std::vector<bool> &)> &,
           const unsigned int);

    template bool
    all_reduce(const bool &,
               const MPI_Comm &,
               const std::function<bool(const bool &, const bool &)> &);

    template std::vector<bool>
    all_reduce(
      const std::vector<bool> &,
      const MPI_Comm &,
      const std::function<std::vector<bool>(const std::vector<bool> &,
                                            const std::vector<bool> &)> &);

    // We need an explicit instantiation of this for the same reason as the
    // other types described in mpi.inst.in
    template void
    internal::all_reduce<bool>(const MPI_Op &,
                               const ArrayView<const bool> &,
                               const MPI_Comm &,
                               const ArrayView<bool> &);


    template bool
    logical_or<bool>(const bool &, const MPI_Comm &);


    template void
    logical_or<bool>(const ArrayView<const bool> &,
                     const MPI_Comm &,
                     const ArrayView<bool> &);


    template std::vector<unsigned int>
    compute_set_union(const std::vector<unsigned int> &vec,
                      const MPI_Comm &                 comm);


    template std::set<unsigned int>
    compute_set_union(const std::set<unsigned int> &set, const MPI_Comm &comm);
#endif

#include "mpi.inst"
  } // end of namespace MPI
} // end of namespace Utilities

DEAL_II_NAMESPACE_CLOSE
