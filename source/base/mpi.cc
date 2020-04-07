// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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
#include <deal.II/base/mpi_tags.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include <iostream>
#include <numeric>
#include <set>
#include <vector>

#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_WITH_MPI
#    include <deal.II/lac/trilinos_parallel_block_vector.h>
#    include <deal.II/lac/trilinos_vector.h>
#    include <deal.II/lac/vector_memory.h>

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
  namespace MPI
  {
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
#  if DEAL_II_MPI_VERSION_GTE(3, 0)
      return MPI_Comm_create_group(comm, group, tag, new_comm);
#  else
      int rank;
      int ierr = MPI_Comm_rank(comm, &rank);
      AssertThrowMPI(ierr);

      int grp_rank;
      ierr = MPI_Group_rank(group, &grp_rank);
      AssertThrowMPI(ierr);
      if (grp_rank == MPI_UNDEFINED)
        {
          *new_comm = MPI_COMM_NULL;
          return MPI_SUCCESS;
        }

      int grp_size;
      ierr = MPI_Group_size(group, &grp_size);
      AssertThrowMPI(ierr);

      ierr = MPI_Comm_dup(MPI_COMM_SELF, new_comm);
      AssertThrowMPI(ierr);

      MPI_Group parent_grp;
      ierr = MPI_Comm_group(comm, &parent_grp);
      AssertThrowMPI(ierr);

      std::vector<int> pids(grp_size);
      std::vector<int> grp_pids(grp_size);
      std::iota(grp_pids.begin(), grp_pids.end(), 0);
      ierr = MPI_Group_translate_ranks(
        group, grp_size, grp_pids.data(), parent_grp, pids.data());
      AssertThrowMPI(ierr);
      ierr = MPI_Group_free(&parent_grp);
      AssertThrowMPI(ierr);

      MPI_Comm comm_old = *new_comm;
      MPI_Comm ic;
      for (int merge_sz = 1; merge_sz < grp_size; merge_sz *= 2)
        {
          const int gid = grp_rank / merge_sz;
          comm_old      = *new_comm;
          if (gid % 2 == 0)
            {
              if ((gid + 1) * merge_sz < grp_size)
                {
                  ierr = (MPI_Intercomm_create(
                    *new_comm, 0, comm, pids[(gid + 1) * merge_sz], tag, &ic));
                  AssertThrowMPI(ierr);
                  ierr = MPI_Intercomm_merge(ic, 0 /* LOW */, new_comm);
                  AssertThrowMPI(ierr);
                }
            }
          else
            {
              ierr = MPI_Intercomm_create(
                *new_comm, 0, comm, pids[(gid - 1) * merge_sz], tag, &ic);
              AssertThrowMPI(ierr);
              ierr = MPI_Intercomm_merge(ic, 1 /* HIGH */, new_comm);
              AssertThrowMPI(ierr);
            }
          if (*new_comm != comm_old)
            {
              ierr = MPI_Comm_free(&ic);
              AssertThrowMPI(ierr);
              ierr = MPI_Comm_free(&comm_old);
              AssertThrowMPI(ierr);
            }
        }

      return MPI_SUCCESS;
#  endif
    }



    std::vector<IndexSet>
    create_ascending_partitioning(const MPI_Comm &           comm,
                                  const IndexSet::size_type &local_size)
    {
      const unsigned int                     n_proc = n_mpi_processes(comm);
      const std::vector<IndexSet::size_type> sizes =
        all_gather(comm, local_size);
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



    /**
     * A re-implementation of compute_point_to_point_communication_pattern
     * using a ConsensusAlgorithm.
     */
    class ConsensusAlgorithmsProcessTargets
      : public ConsensusAlgorithms::Process<unsigned int, unsigned int>
    {
    public:
      ConsensusAlgorithmsProcessTargets(const std::vector<unsigned int> &target)
        : target(target)
      {}

      using T1 = unsigned int;
      using T2 = unsigned int;

      virtual void
      answer_request(const unsigned int other_rank,
                     const std::vector<T1> &,
                     std::vector<T2> &) override
      {
        this->sources.push_back(other_rank);
      }

      /**
       * Simply return the user-provided list.
       *
       * @return List of processes this process wants to send requests to.
       */
      virtual std::vector<unsigned int>
      compute_targets() override
      {
        return target;
      }

      /**
       * The result of the consensus algorithm.
       * @return Sorted list of ranks of processes wanting to send a request to
       *         this process.
       */
      std::vector<unsigned int>
      get_result()
      {
        std::sort(sources.begin(), sources.end());
        return sources;
      }

    private:
      /**
       * List of processes this process wants to send requests to.
       */
      const std::vector<unsigned int> &target;

      /**
       * List of ranks of processes wanting to send a request to this process.
       */
      std::vector<unsigned int> sources;
    };



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
          Assert(destination != myid,
                 ExcMessage(
                   "There is no point in communicating with ourselves."));
        }

#  if DEAL_II_MPI_VERSION_GTE(3, 0)

      ConsensusAlgorithmsProcessTargets process(destinations);
      ConsensusAlgorithms::NBX<ConsensusAlgorithmsProcessTargets::T1,
                               ConsensusAlgorithmsProcessTargets::T2>
        consensus_algorithm(process, mpi_comm);
      consensus_algorithm.run();
      return process.get_result();

#  elif DEAL_II_MPI_VERSION_GTE(2, 2)

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
#  else
      // let all processors communicate the maximal number of destinations
      // they have
      const unsigned int max_n_destinations =
        Utilities::MPI::max(destinations.size(), mpi_comm);

      if (max_n_destinations == 0)
        // all processes have nothing to send/receive:
        return std::vector<unsigned int>();

      // now that we know the number of data packets every processor wants to
      // send, set up a buffer with the maximal size and copy our destinations
      // in there, padded with -1's
      std::vector<unsigned int> my_destinations(max_n_destinations,
                                                numbers::invalid_unsigned_int);
      std::copy(destinations.begin(),
                destinations.end(),
                my_destinations.begin());

      // now exchange these (we could communicate less data if we used
      // MPI_Allgatherv, but we'd have to communicate my_n_destinations to all
      // processors in this case, which is more expensive than the reduction
      // operation above in MPI_Allreduce)
      std::vector<unsigned int> all_destinations(max_n_destinations * n_procs);
      const int                 ierr = MPI_Allgather(my_destinations.data(),
                                     max_n_destinations,
                                     MPI_UNSIGNED,
                                     all_destinations.data(),
                                     max_n_destinations,
                                     MPI_UNSIGNED,
                                     mpi_comm);
      AssertThrowMPI(ierr);

      // now we know who is going to communicate with whom. collect who is
      // going to communicate with us!
      std::vector<unsigned int> origins;
      for (unsigned int i = 0; i < n_procs; ++i)
        for (unsigned int j = 0; j < max_n_destinations; ++j)
          if (all_destinations[i * max_n_destinations + j] == myid)
            origins.push_back(i);
          else if (all_destinations[i * max_n_destinations + j] ==
                   numbers::invalid_unsigned_int)
            break;

      return origins;
#  endif
    }



    unsigned int
    compute_n_point_to_point_communications(
      const MPI_Comm &                 mpi_comm,
      const std::vector<unsigned int> &destinations)
    {
      const unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_comm);

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

#  if DEAL_II_MPI_VERSION_GTE(2, 2)
      // Find out how many processes will send to this one
      // MPI_Reduce_scatter(_block) does exactly this
      unsigned int n_recv_from = 0;

      const int ierr = MPI_Reduce_scatter_block(
        dest_vector.data(), &n_recv_from, 1, MPI_UNSIGNED, MPI_SUM, mpi_comm);

      AssertThrowMPI(ierr);

      return n_recv_from;
#  else
      // Find out how many processes will send to this one
      // by reducing with sum and then scattering the
      // results over all processes
      std::vector<unsigned int> buffer(dest_vector.size());
      unsigned int              n_recv_from = 0;

      MPI_Reduce(dest_vector.data(),
                 buffer.data(),
                 dest_vector.size(),
                 MPI_UNSIGNED,
                 MPI_SUM,
                 0,
                 mpi_comm);
      MPI_Scatter(buffer.data(),
                  1,
                  MPI_UNSIGNED,
                  &n_recv_from,
                  1,
                  MPI_UNSIGNED,
                  0,
                  mpi_comm);

      return n_recv_from;
#  endif
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

        for (int i = 0; i < *len; i++)
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
          for (unsigned int i = 0; i < my_values.size(); i++)
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

      AssertDimension(Utilities::MPI::min(my_values.size(), mpi_communicator),
                      Utilities::MPI::max(my_values.size(), mpi_communicator));

      AssertDimension(my_values.size(), result.size());



      // To avoid uninitialized values on some MPI implementations, provide
      // result with a default value already...
      MinMaxAvg dummy = {0.,
                         std::numeric_limits<double>::max(),
                         -std::numeric_limits<double>::max(),
                         0,
                         0,
                         0.};

      for (auto &i : result)
        i = dummy;

      const unsigned int my_id =
        dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
      const unsigned int numproc =
        dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);

      MPI_Op op;
      int    ierr =
        MPI_Op_create(reinterpret_cast<MPI_User_function *>(&max_reduce),
                      true,
                      &op);
      AssertThrowMPI(ierr);

      std::vector<MinMaxAvg> in(my_values.size());

      for (unsigned int i = 0; i < my_values.size(); i++)
        {
          in[i].sum = in[i].min = in[i].max = my_values[i];
          in[i].min_index = in[i].max_index = my_id;
        }

      MPI_Datatype type;
      int          lengths[]       = {3, 2, 1};
      MPI_Aint     displacements[] = {0,
                                  offsetof(MinMaxAvg, min_index),
                                  offsetof(MinMaxAvg, avg)};
      MPI_Datatype types[]         = {MPI_DOUBLE, MPI_INT, MPI_DOUBLE};

      ierr = MPI_Type_create_struct(3, lengths, displacements, types, &type);
      AssertThrowMPI(ierr);

      ierr = MPI_Type_commit(&type);
      AssertThrowMPI(ierr);
      ierr = MPI_Allreduce(
        in.data(), result.data(), my_values.size(), type, op, mpi_communicator);
      AssertThrowMPI(ierr);

      ierr = MPI_Type_free(&type);
      AssertThrowMPI(ierr);

      ierr = MPI_Op_free(&op);
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



    std::vector<IndexSet>
    create_ascending_partitioning(const MPI_Comm & /*comm*/,
                                  const IndexSet::size_type &local_size)
    {
      return std::vector<IndexSet>(1, complete_index_set(local_size));
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

      for (unsigned int i = 0; i < my_values.size(); i++)
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
      ierr = SlepcInitialize(&argc, &argv, nullptr, nullptr);
      AssertThrow(ierr == 0, SLEPcWrappers::SolverBase::ExcSLEPcError(ierr));
#  else
      // or just initialize PETSc alone:
      ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
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
#  if DEAL_II_P4EST_VERSION_GTE(2, 0, 0, 0)
#  else
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
#  if defined(DEAL_II_WITH_TRILINOS)
      GrowingVectorMemory<
        TrilinosWrappers::MPI::Vector>::release_unused_memory();
      GrowingVectorMemory<
        TrilinosWrappers::MPI::BlockVector>::release_unused_memory();
#  endif
#endif


      // Now deal with PETSc (with or without MPI). Only delete the vectors if
      // finalize hasn't been called yet, otherwise this will lead to errors.
#ifdef DEAL_II_WITH_PETSC
      if ((PetscInitializeCalled == PETSC_TRUE) &&
          (PetscFinalizeCalled == PETSC_FALSE))
        {
          GrowingVectorMemory<
            PETScWrappers::MPI::Vector>::release_unused_memory();
          GrowingVectorMemory<
            PETScWrappers::MPI::BlockVector>::release_unused_memory();

#  ifdef DEAL_II_WITH_SLEPC
          // and now end SLEPc (with PETSc)
          SlepcFinalize();
#  else
          // or just end PETSc.
          PetscFinalize();
#  endif
        }
#endif

// There is a similar issue with CUDA: The destructor of static objects might
// run after the CUDA driver is unloaded. Hence, also release all memory
// related to CUDA vectors.
#ifdef DEAL_II_WITH_CUDA
      GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>>::
        release_unused_memory();
      GrowingVectorMemory<
        LinearAlgebra::distributed::Vector<float, MemorySpace::CUDA>>::
        release_unused_memory();
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
        std::pair<types::global_dof_index, types::global_dof_index>,
        unsigned int>
        consensus_algorithm(process, comm);
      consensus_algorithm.run();

      return owning_ranks;
    }



    CollectiveMutex::CollectiveMutex()
      : locked(false)
      , request(MPI_REQUEST_NULL)
    {
      Utilities::MPI::MPI_InitFinalize::register_request(request);
    }



    CollectiveMutex::~CollectiveMutex()
    {
      Assert(
        !locked,
        ExcMessage(
          "Error: MPI::CollectiveMutex is still locked while being destroyed!"));

      Utilities::MPI::MPI_InitFinalize::unregister_request(request);
    }



    void
    CollectiveMutex::lock(MPI_Comm comm)
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

#  if 0 && DEAL_II_MPI_VERSION_GTE(3, 0)
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
    CollectiveMutex::unlock(MPI_Comm comm)
    {
      (void)comm;

      Assert(
        locked,
        ExcMessage(
          "Error: MPI::CollectiveMutex needs to be locked before unlock()"));

#ifdef DEAL_II_WITH_MPI

      // TODO: For now, we implement this mutex with a blocking barrier
      // in the lock and unlock. It needs to be tested, if we can move
      // to a nonblocking barrier (code disabled below):

#  if 0 && DEAL_II_MPI_VERSION_GTE(3, 0)
      const int ierr = MPI_Ibarrier(comm, &request);
      AssertThrowMPI(ierr);
#  else
      const int ierr = MPI_Barrier(comm);
      AssertThrowMPI(ierr);
#  endif
#endif

      locked = false;
    }


    template std::vector<unsigned int>
    compute_set_union(const std::vector<unsigned int> &vec,
                      const MPI_Comm &                 comm);


    template std::set<unsigned int>
    compute_set_union(const std::set<unsigned int> &set, const MPI_Comm &comm);

#include "mpi.inst"
  } // end of namespace MPI
} // end of namespace Utilities

DEAL_II_NAMESPACE_CLOSE
