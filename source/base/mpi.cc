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



    std::vector<unsigned int>
    compute_point_to_point_communication_pattern(
      const MPI_Comm &                 mpi_comm,
      const std::vector<unsigned int> &destinations)
    {
      const unsigned int myid    = Utilities::MPI::this_mpi_process(mpi_comm);
      const unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_comm);

      for (const unsigned int destination : destinations)
        {
          (void)destination;
          Assert(destination < n_procs, ExcIndexRange(destination, 0, n_procs));
          Assert(destination != myid,
                 ExcMessage(
                   "There is no point in communicating with ourselves."));
        }

#  if DEAL_II_MPI_VERSION_GTE(2, 2)
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
        MPI_Isend(&myid,
                  1,
                  MPI_UNSIGNED,
                  el,
                  32766,
                  mpi_comm,
                  send_requests.data() + (&el - destinations.data()));

      // if no one to receive from, return an empty vector
      if (n_recv_from == 0)
        return std::vector<unsigned int>();

      // ...otherwise receive `n_recv_from` times from the processes
      // who communicate with this one. Store the obtained id's
      // in the resulting vector
      std::vector<unsigned int> origins(n_recv_from);
      for (auto &el : origins)
        MPI_Recv(&el,
                 1,
                 MPI_UNSIGNED,
                 MPI_ANY_SOURCE,
                 32766,
                 mpi_comm,
                 MPI_STATUS_IGNORE);

      MPI_Waitall(destinations.size(),
                  send_requests.data(),
                  MPI_STATUSES_IGNORE);
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
          Assert(destination < n_procs, ExcIndexRange(destination, 0, n_procs));
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
        (void)len;
        const MinMaxAvg *in_lhs    = static_cast<const MinMaxAvg *>(in_lhs_);
        MinMaxAvg *      inout_rhs = static_cast<MinMaxAvg *>(inout_rhs_);

        Assert(*len == 1, ExcInternalError());

        inout_rhs->sum += in_lhs->sum;
        if (inout_rhs->min > in_lhs->min)
          {
            inout_rhs->min       = in_lhs->min;
            inout_rhs->min_index = in_lhs->min_index;
          }
        else if (inout_rhs->min == in_lhs->min)
          {
            // choose lower cpu index when tied to make operator commutative
            if (inout_rhs->min_index > in_lhs->min_index)
              inout_rhs->min_index = in_lhs->min_index;
          }

        if (inout_rhs->max < in_lhs->max)
          {
            inout_rhs->max       = in_lhs->max;
            inout_rhs->max_index = in_lhs->max_index;
          }
        else if (inout_rhs->max == in_lhs->max)
          {
            // choose lower cpu index when tied to make operator commutative
            if (inout_rhs->max_index > in_lhs->max_index)
              inout_rhs->max_index = in_lhs->max_index;
          }
      }
    } // namespace



    MinMaxAvg
    min_max_avg(const double my_value, const MPI_Comm &mpi_communicator)
    {
      // If MPI was not started, we have a serial computation and cannot run
      // the other MPI commands
      if (job_supports_mpi() == false)
        {
          MinMaxAvg result;
          result.sum       = my_value;
          result.avg       = my_value;
          result.min       = my_value;
          result.max       = my_value;
          result.min_index = 0;
          result.max_index = 0;

          return result;
        }

      // To avoid uninitialized values on some MPI implementations, provide
      // result with a default value already...
      MinMaxAvg result = {0.,
                          std::numeric_limits<double>::max(),
                          -std::numeric_limits<double>::max(),
                          0,
                          0,
                          0.};

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

      MinMaxAvg in;
      in.sum = in.min = in.max = my_value;
      in.min_index = in.max_index = my_id;

      MPI_Datatype type;
      int          lengths[]       = {3, 2};
      MPI_Aint     displacements[] = {0, offsetof(MinMaxAvg, min_index)};
      MPI_Datatype types[]         = {MPI_DOUBLE, MPI_INT};

      ierr = MPI_Type_create_struct(2, lengths, displacements, types, &type);
      AssertThrowMPI(ierr);

      ierr = MPI_Type_commit(&type);
      AssertThrowMPI(ierr);
      ierr = MPI_Allreduce(&in, &result, 1, type, op, mpi_communicator);
      AssertThrowMPI(ierr);

      ierr = MPI_Type_free(&type);
      AssertThrowMPI(ierr);

      ierr = MPI_Op_free(&op);
      AssertThrowMPI(ierr);

      result.avg = result.sum / numproc;

      return result;
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



    MinMaxAvg
    min_max_avg(const double my_value, const MPI_Comm &)
    {
      MinMaxAvg result;

      result.sum       = my_value;
      result.avg       = my_value;
      result.min       = my_value;
      result.max       = my_value;
      result.min_index = 0;
      result.max_index = 0;

      return result;
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


    MPI_InitFinalize::~MPI_InitFinalize()
    {
      // make memory pool release all PETSc/Trilinos/MPI-based vectors that
      // are no longer used at this point. this is relevant because the static
      // object destructors run for these vectors at the end of the program
      // would run after MPI_Finalize is called, leading to errors

#ifdef DEAL_II_WITH_MPI
      // Start with the deal.II MPI vectors (need to do this before finalizing
      // PETSc because it finalizes MPI).  Delete vectors from the pools:
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
              std::cerr
                << "ERROR: Uncaught exception in MPI_InitFinalize on proc "
                << this_mpi_process(MPI_COMM_WORLD)
                << ". Skipping MPI_Finalize() to avoid a deadlock."
                << std::endl;
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

    template <typename T1, typename T2>
    void
    ConsensusAlgorithmProcess<T1, T2>::process_request(const unsigned int,
                                                       const std::vector<T1> &,
                                                       std::vector<T2> &)
    {
      // noting to do
    }



    template <typename T1, typename T2>
    void
    ConsensusAlgorithmProcess<T1, T2>::pack_recv_buffer(const int,
                                                        std::vector<T1> &)
    {
      // noting to do
    }



    template <typename T1, typename T2>
    void
    ConsensusAlgorithmProcess<T1, T2>::prepare_recv_buffer(const int,
                                                           std::vector<T2> &)
    {
      // noting to do
    }



    template <typename T1, typename T2>
    void
    ConsensusAlgorithmProcess<T1, T2>::unpack_recv_buffer(
      const int,
      const std::vector<T2> &)
    {
      // noting to do
    }



    template <typename T1, typename T2>
    ConsensusAlgorithm<T1, T2>::ConsensusAlgorithm(
      ConsensusAlgorithmProcess<T1, T2> &process,
      const MPI_Comm &                   comm)
      : process(process)
      , comm(comm)
      , my_rank(this_mpi_process(comm))
      , n_procs(n_mpi_processes(comm))
    {}



    template <typename T1, typename T2>
    ConsensusAlgorithm_NBX<T1, T2>::ConsensusAlgorithm_NBX(
      ConsensusAlgorithmProcess<T1, T2> &process,
      const MPI_Comm &                   comm)
      : ConsensusAlgorithm<T1, T2>(process, comm)
    {}



    template <typename T1, typename T2>
    void
    ConsensusAlgorithm_NBX<T1, T2>::run()
    {
      // 1) send requests and start receiving the answers
      start_communication();

      // 2) answer requests and check if all requests of this process have been
      //    answered
      while (!check_own_state())
        process_requests();

      // 3) signal to all other processes that all requests of this process have
      //    been answered
      signal_finish();

      // 4) nevertheless, this process has to keep on answering (potential)
      //    incoming requests until all processes have received the
      //    answer to all requests
      while (!check_global_state())
        process_requests();

      // 5) process the answer to all requests
      clean_up_and_end_communication();
    }



    template <typename T1, typename T2>
    bool
    ConsensusAlgorithm_NBX<T1, T2>::check_own_state()
    {
#ifdef DEAL_II_WITH_MPI
      int        all_receive_requests_are_done;
      const auto ierr = MPI_Testall(recv_requests.size(),
                                    recv_requests.data(),
                                    &all_receive_requests_are_done,
                                    MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);

      return all_receive_requests_are_done;
#else
      return true;
#endif
    }



    template <typename T1, typename T2>
    void
    ConsensusAlgorithm_NBX<T1, T2>::signal_finish()
    {
#ifdef DEAL_II_WITH_MPI
#  if DEAL_II_MPI_VERSION_GTE(3, 0)
      const auto ierr = MPI_Ibarrier(this->comm, &barrier_request);
      AssertThrowMPI(ierr);
#  else
      AssertThrow(
        false,
        ExcMessage(
          "ConsensusAlgorithm_NBX uses MPI 3.0 features. You should compile with at least MPI 3.0."));
#  endif
#endif
    }



    template <typename T1, typename T2>
    bool
    ConsensusAlgorithm_NBX<T1, T2>::check_global_state()
    {
#ifdef DEAL_II_WITH_MPI
      int        all_ranks_reached_barrier;
      const auto ierr = MPI_Test(&barrier_request,
                                 &all_ranks_reached_barrier,
                                 MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
      return all_ranks_reached_barrier;
#else
      return true;
#endif
    }



    template <typename T1, typename T2>
    void
    ConsensusAlgorithm_NBX<T1, T2>::process_requests()
    {
#ifdef DEAL_II_WITH_MPI
      // check if there is a request pending
      MPI_Status status;
      int        request_is_pending;
      const auto ierr = MPI_Iprobe(
        MPI_ANY_SOURCE, tag_request, this->comm, &request_is_pending, &status);
      AssertThrowMPI(ierr);

      if (request_is_pending) // request is pending
        {
          // get rank of requesting process
          const auto other_rank = status.MPI_SOURCE;

#  ifdef DEBUG
          Assert(requesting_processes.find(other_rank) ==
                   requesting_processes.end(),
                 ExcMessage("Process is requesting a second time!"));
          requesting_processes.insert(other_rank);
#  endif

          std::vector<T1> buffer_recv;
          // get size of of incoming message
          int  number_amount;
          auto ierr = MPI_Get_count(&status, MPI_BYTE, &number_amount);
          AssertThrowMPI(ierr);

          // allocate memory for incoming message
          Assert(number_amount % sizeof(T1) == 0, ExcInternalError());
          buffer_recv.resize(number_amount / sizeof(T1));
          ierr = MPI_Recv(buffer_recv.data(),
                          number_amount,
                          MPI_BYTE,
                          other_rank,
                          tag_request,
                          this->comm,
                          &status);
          AssertThrowMPI(ierr);

          // allocate memory for answer message
          request_buffers.emplace_back();
          request_requests.emplace_back(new MPI_Request);

          // process request
          auto &request_buffer = request_buffers.back();
          this->process.process_request(other_rank,
                                        buffer_recv,
                                        request_buffer);

          // start to send answer back
          ierr = MPI_Isend(request_buffer.data(),
                           request_buffer.size() * sizeof(T2),
                           MPI_BYTE,
                           other_rank,
                           tag_delivery,
                           this->comm,
                           request_requests.back().get());
          AssertThrowMPI(ierr);
        }
#endif
    }



    template <typename T1, typename T2>
    void
    ConsensusAlgorithm_NBX<T1, T2>::start_communication()
    {
#ifdef DEAL_II_WITH_MPI
      // 1)
      targets              = this->process.compute_targets();
      const auto n_targets = targets.size();

      // 2) allocate memory
      recv_buffers.resize(n_targets);
      recv_requests.resize(n_targets);
      send_requests.resize(n_targets);
      send_buffers.resize(n_targets);

      {
        // 4) send and receive
        for (unsigned int i = 0; i < n_targets; i++)
          {
            const unsigned int rank  = targets[i];
            const unsigned int index = i;

            // translate index set to a list of pairs
            auto &send_buffer = send_buffers[index];
            this->process.pack_recv_buffer(rank, send_buffer);

            // start to send data
            auto ierr = MPI_Isend(send_buffer.data(),
                                  send_buffer.size() * sizeof(T1),
                                  MPI_BYTE,
                                  rank,
                                  tag_request,
                                  this->comm,
                                  &send_requests[index]);
            AssertThrowMPI(ierr);

            // start to receive data
            auto &recv_buffer = recv_buffers[index];
            this->process.prepare_recv_buffer(rank, recv_buffer);
            ierr = MPI_Irecv(recv_buffer.data(),
                             recv_buffer.size() * sizeof(T2),
                             MPI_BYTE,
                             rank,
                             tag_delivery,
                             this->comm,
                             &recv_requests[index]);
            AssertThrowMPI(ierr);
          }
      }
#endif
    }



    template <typename T1, typename T2>
    void
    ConsensusAlgorithm_NBX<T1, T2>::clean_up_and_end_communication()
    {
#ifdef DEAL_II_WITH_MPI
      // clean up
      {
        auto ierr = MPI_Waitall(send_requests.size(),
                                send_requests.data(),
                                MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr);

        ierr = MPI_Waitall(recv_requests.size(),
                           recv_requests.data(),
                           MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr);

        ierr = MPI_Wait(&barrier_request, MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);

        for (auto &i : request_requests)
          {
            const auto ierr = MPI_Wait(i.get(), MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);
          }

#  ifdef DEBUG
        // note: IBarrier seems to make problem during testing, this additional
        // Barrier seems to help
        MPI_Barrier(this->comm);
#  endif
      }

      // unpack data
      {
        for (unsigned int i = 0; i < targets.size(); i++)
          this->process.unpack_recv_buffer(targets[i], recv_buffers[i]);
      }
#endif
    }



    /**
     * A re-implementation of compute_point_to_point_communication_pattern
     * using the ConsensusAlgorithm.
     */
    class ConsensusAlgorithmProcessTargets
      : public ConsensusAlgorithmProcess<int, int>
    {
    public:
      ConsensusAlgorithmProcessTargets(std::vector<unsigned int> &target)
        : target(target)
      {}

      using T1 = int;
      using T2 = int;

      virtual void
      process_request(const unsigned int other_rank,
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



    template <typename T1, typename T2>
    ConsensusAlgorithm_PEX<T1, T2>::ConsensusAlgorithm_PEX(
      ConsensusAlgorithmProcess<T1, T2> &process,
      const MPI_Comm &                   comm)
      : ConsensusAlgorithm<T1, T2>(process, comm)
    {}



    template <typename T1, typename T2>
    void
    ConsensusAlgorithm_PEX<T1, T2>::run()
    {
      // 1) send requests and start receiving the answers
      //    especially determine how many requests are expected
      const unsigned int n_requests = start_communication();

      // 2) answer requests
      for (unsigned int request = 0; request < n_requests; request++)
        process_requests(request);

      // 3) process answers
      clean_up_and_end_communication();
    }



    template <typename T1, typename T2>
    void
    ConsensusAlgorithm_PEX<T1, T2>::process_requests(int index)
    {
#ifdef DEAL_II_WITH_MPI
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, tag_request, this->comm, &status);

      // get rank of incoming message
      const auto other_rank = status.MPI_SOURCE;

      std::vector<T1> buffer_recv;

      // get size of incoming message
      int  number_amount;
      auto ierr = MPI_Get_count(&status, MPI_BYTE, &number_amount);
      AssertThrowMPI(ierr);

      // allocate memory for incoming message
      Assert(number_amount % sizeof(T1) == 0, ExcInternalError());
      buffer_recv.resize(number_amount / sizeof(T1));
      ierr = MPI_Recv(buffer_recv.data(),
                      number_amount,
                      MPI_BYTE,
                      other_rank,
                      tag_request,
                      this->comm,
                      &status);
      AssertThrowMPI(ierr);

      // process request
      auto &request_buffer = requests_buffers[index];
      this->process.process_request(other_rank, buffer_recv, request_buffer);

      // start to send answer back
      ierr = MPI_Isend(request_buffer.data(),
                       request_buffer.size() * sizeof(T2),
                       MPI_BYTE,
                       other_rank,
                       tag_delivery,
                       this->comm,
                       &requests_answers[index]);
      AssertThrowMPI(ierr);
#else
      (void)index;
#endif
    }



    template <typename T1, typename T2>
    unsigned int
    ConsensusAlgorithm_PEX<T1, T2>::start_communication()
    {
#ifdef DEAL_II_WITH_MPI
      // 1) determine with which processes this process wants to communicate
      targets = this->process.compute_targets();

      // 2) determine who wants to communicate with this process
      const bool use_nbx = false;
      if (!use_nbx)
        {
          sources =
            compute_point_to_point_communication_pattern(this->comm, targets);
        }
      else
        {
          ConsensusAlgorithmProcessTargets process(targets);
          ConsensusAlgorithm_NBX<ConsensusAlgorithmProcessTargets::T1,
                                 ConsensusAlgorithmProcessTargets::T2>
            consensus_algorithm(process, this->comm);
          consensus_algorithm.run();
          sources = process.get_result();
        }

      const auto n_targets = targets.size();
      const auto n_sources = sources.size();

      // 2) allocate memory
      recv_buffers.resize(n_targets);
      send_buffers.resize(n_targets);
      send_and_recv_buffers.resize(2 * n_targets);

      requests_answers.resize(n_sources);
      requests_buffers.resize(n_sources);

      // 4) send and receive
      for (unsigned int i = 0; i < n_targets; i++)
        {
          const unsigned int rank = targets[i];

          // pack data which should be sent
          auto &send_buffer = send_buffers[i];
          this->process.pack_recv_buffer(rank, send_buffer);

          // start to send data
          auto ierr = MPI_Isend(send_buffer.data(),
                                send_buffer.size() * sizeof(T1),
                                MPI_BYTE,
                                rank,
                                tag_request,
                                this->comm,
                                &send_and_recv_buffers[n_targets + i]);
          AssertThrowMPI(ierr);

          // start to receive data
          auto &recv_buffer = recv_buffers[i];
          this->process.prepare_recv_buffer(rank, recv_buffer);
          ierr = MPI_Irecv(recv_buffer.data(),
                           recv_buffer.size() * sizeof(T2),
                           MPI_BYTE,
                           rank,
                           tag_delivery,
                           this->comm,
                           &send_and_recv_buffers[i]);
          AssertThrowMPI(ierr);
        }

      return sources.size();
#else
      return 0;
#endif
    }



    template <typename T1, typename T2>
    void
    ConsensusAlgorithm_PEX<T1, T2>::clean_up_and_end_communication()
    {
#ifdef DEAL_II_WITH_MPI
      // finalize all MPI_Requests
      MPI_Waitall(send_and_recv_buffers.size(),
                  send_and_recv_buffers.data(),
                  MPI_STATUSES_IGNORE);
      MPI_Waitall(requests_answers.size(),
                  requests_answers.data(),
                  MPI_STATUSES_IGNORE);

      // unpack received data
      for (unsigned int i = 0; i < targets.size(); i++)
        this->process.unpack_recv_buffer(targets[i], recv_buffers[i]);
#endif
    }



    template <typename T1, typename T2>
    ConsensusAlgorithmSelector<T1, T2>::ConsensusAlgorithmSelector(
      ConsensusAlgorithmProcess<T1, T2> &process,
      const MPI_Comm &                   comm)
      : ConsensusAlgorithm<T1, T2>(process, comm)
    {
      // Depending on the number of processes we switch between implementations.
      // We reduce the threshold for debug mode to be able to test also the
      // non-blocking implementation. This feature is tested by:
      // tests/multigrid/transfer_matrix_free_06.with_mpi=true.with_p4est=true.with_trilinos=true.mpirun=15.output
#ifdef DEAL_II_WITH_MPI
#  if DEAL_II_MPI_VERSION_GTE(3, 0)
#    ifdef DEBUG
      if (Utilities::MPI::n_mpi_processes(comm) > 14)
#    else
      if (Utilities::MPI::n_mpi_processes(comm) > 99)
#    endif
        consensus_algo.reset(new ConsensusAlgorithm_NBX<T1, T2>(process, comm));
      else
#  endif
#endif
        consensus_algo.reset(new ConsensusAlgorithm_PEX<T1, T2>(process, comm));
    }



    template <typename T1, typename T2>
    void
    ConsensusAlgorithmSelector<T1, T2>::run()
    {
      consensus_algo->run();
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
      internal::ComputeIndexOwner::ConsensusAlgorithmPayload process(
        owned_indices, indices_to_look_up, comm, owning_ranks);

      // Step 2: read dictionary
      // Communicate with the process who owns the index in the static
      // partition (i.e. in the dictionary). This process returns the actual
      // owner of the index.
      ConsensusAlgorithmSelector<
        std::pair<types::global_dof_index, types::global_dof_index>,
        unsigned int>
        consensus_algorithm(process, comm);
      consensus_algorithm.run();

      return owning_ranks;
    }

    template class ConsensusAlgorithmSelector<
      std::pair<types::global_dof_index, types::global_dof_index>,
      unsigned int>;

#include "mpi.inst"
  } // end of namespace MPI
} // end of namespace Utilities

DEAL_II_NAMESPACE_CLOSE
