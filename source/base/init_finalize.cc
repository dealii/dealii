// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2025 by the deal.II authors
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

#include <deal.II/base/init_finalize.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>

#include <Kokkos_Core.hpp>

#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_WITH_MPI
#    include <deal.II/lac/trilinos_parallel_block_vector.h>
#    include <deal.II/lac/trilinos_vector.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#    include <Epetra_MpiComm.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
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

#include <set>
#include <string>


DEAL_II_NAMESPACE_OPEN


/* Force initialization of static struct: */
InitFinalize::Signals InitFinalize::signals = InitFinalize::Signals();


InitFinalize::InitFinalize([[maybe_unused]] int    &argc,
                           [[maybe_unused]] char **&argv,
                           const InitializeLibrary &libraries,
                           const unsigned int       max_num_threads)
  : libraries(libraries)
{
  [[maybe_unused]] static bool constructor_has_already_run = false;
  Assert(constructor_has_already_run == false,
         ExcMessage("You can only create a single object of this class "
                    "in a program since it initializes the MPI system."));


  [[maybe_unused]] int ierr = 0;
#ifdef DEAL_II_WITH_MPI
  if (static_cast<bool>(libraries & InitializeLibrary::MPI))
    {
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
    }
#endif

    // we are allowed to call MPI_Init ourselves and PETScInitialize will
    // detect this. This allows us to use MPI_Init_thread instead.
#ifdef DEAL_II_WITH_PETSC
  PetscErrorCode pierr;
#  ifdef DEAL_II_WITH_SLEPC
  // Initialize SLEPc (with PETSc):
  if (static_cast<bool>(libraries & InitializeLibrary::SLEPc))
    {
      finalize_petscslepc = SlepcInitializeCalled ? false : true;
      pierr               = SlepcInitialize(&argc, &argv, nullptr, nullptr);
      AssertThrow(pierr == 0, SLEPcWrappers::SolverBase::ExcSLEPcError(pierr));
    }
#  else
  // or just initialize PETSc alone:
  if (static_cast<bool>(libraries & InitializeLibrary::PETSc))
    {
      finalize_petscslepc = PetscInitializeCalled ? false : true;
      pierr               = PetscInitialize(&argc, &argv, nullptr, nullptr);
      AssertThrow(pierr == 0, ExcPETScError(pierr));
    }
#  endif

  // Disable PETSc exception handling. This just prints a large wall
  // of text that is not particularly helpful for what we do:
  if (static_cast<bool>(libraries & InitializeLibrary::SLEPc) ||
      static_cast<bool>(libraries & InitializeLibrary::PETSc))
    {
      pierr = PetscPopSignalHandler();
      AssertThrow(pierr == 0, ExcPETScError(pierr));
    }
#endif

    // Initialize zoltan
#ifdef DEAL_II_TRILINOS_WITH_ZOLTAN
  if (static_cast<bool>(libraries & InitializeLibrary::Zoltan))
    {
      float version;
      Zoltan_Initialize(argc, argv, &version);
    }
#endif

    // Initialize p4est and libsc components
#ifdef DEAL_II_WITH_P4EST
  if (static_cast<bool>(libraries & InitializeLibrary::P4EST))
    {
#  if DEAL_II_P4EST_VERSION_GTE(2, 5, 0, 0)
      // This feature is broken in version 2.0.0 for calls to
      // MPI_Comm_create_group (see cburstedde/p4est#30).
      // Disabling it leads to more verbose p4est error messages
      // which should be fine.
      sc_init(MPI_COMM_WORLD, 0, 0, nullptr, SC_LP_SILENT);
#  endif
      p4est_init(nullptr, SC_LP_SILENT);
    }
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
      unsigned int n_threads = MultithreadInfo::n_cores();
#ifdef DEAL_II_WITH_MPI
      if (static_cast<bool>(libraries & InitializeLibrary::MPI))
        {
          int MPI_has_been_started = 0;
          int ierr                 = MPI_Initialized(&MPI_has_been_started);
          AssertThrowMPI(ierr);

          // we need to figure out how many MPI processes there are on the
          // current node, as well as how many CPU cores we have. for the
          // first task, check what get_hostname() returns and then do an
          // allgather so each processor gets the answer
          //
          // in calculating the length of the string, don't forget the
          // terminating \0 on C-style strings
          const std::string hostname = Utilities::System::get_hostname();

          int my_hostname_size  = hostname.size() + 1;
          int max_hostname_size = -1;
          ierr                  = MPI_Allreduce(&my_hostname_size,
                               &max_hostname_size,
                               1,
                               MPI_INT,
                               MPI_MAX,
                               MPI_COMM_WORLD);
          AssertThrowMPI(ierr);
          std::vector<char> hostname_array(max_hostname_size);
          std::copy(hostname.c_str(),
                    hostname.c_str() + hostname.size() + 1,
                    hostname_array.begin());

          int n_mpi_processes = 1;
          if (MPI_has_been_started)
            {
              ierr = MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_processes);
              AssertThrowMPI(ierr);
            }
          std::vector<char> all_hostnames(max_hostname_size * n_mpi_processes);
          ierr = MPI_Allgather(hostname_array.data(),
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
          int          rank                = 0;
          if (MPI_has_been_started)
            {
              ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
              AssertThrowMPI(ierr);
            }
          for (int i = 0; i < n_mpi_processes; ++i)
            if (std::string(all_hostnames.data() + i * max_hostname_size) ==
                hostname)
              {
                ++n_local_processes;
                if (i <= rank)
                  ++nth_process_on_host;
              }
          Assert(nth_process_on_host > 0, ExcInternalError());


          // compute how many cores each process gets. if the number does not
          // divide evenly, then we get one more core if we are among the
          // first few processes
          //
          // if the number would be zero, round up to one since every process
          // needs to have at least one thread
          n_threads =
            std::max(MultithreadInfo::n_cores() / n_local_processes +
                       (nth_process_on_host <=
                            MultithreadInfo::n_cores() % n_local_processes ?
                          1 :
                          0),
                     1U);
        }
#endif

      // finally set this number of threads
      MultithreadInfo::set_thread_limit(n_threads);
    }

  // Initialize Kokkos
  if (static_cast<bool>(libraries & InitializeLibrary::Kokkos))
    {
      // argv has argc+1 elements and the last one is a nullptr. For appending
      // one element we thus create a new argv by copying the first argc
      // elements, append the new option, and then a nullptr.
      //
      // We do get in trouble, though, if a user program is called with
      // '--help' as a command line argument. This '--help' gets passed on to
      // Kokkos, which promptly responds with a lengthy message that the user
      // likely did not intend. As a consequence, filter out this specific
      // flag.
      std::vector<char *> argv_new;
      for (auto *const arg : make_array_view(&argv[0], &argv[0] + argc))
        if (std::strcmp(arg, "--help") != 0)
          argv_new.push_back(arg);

      std::stringstream threads_flag;
#if DEAL_II_KOKKOS_VERSION_GTE(3, 7, 0)
      threads_flag << "--kokkos-num-threads=" << MultithreadInfo::n_threads();
#else
      threads_flag << "--kokkos-threads=" << MultithreadInfo::n_threads();
#endif
      const std::string threads_flag_string = threads_flag.str();
      argv_new.push_back(const_cast<char *>(threads_flag_string.c_str()));
      argv_new.push_back(nullptr);

      // The first argument in Kokkos::initialize is of type int&. Hence, we
      // need to define a new variable to pass to it (instead of using argc+1
      // inline).
      int argc_new = argv_new.size() - 1;
      Kokkos::initialize(argc_new, argv_new.data());
    }

  // As a final step call the at_mpi_init() signal handler.
  signals.at_mpi_init();
}



void
InitFinalize::register_request(MPI_Request &request)
{
  // insert if it is not in the set already:
  requests.insert(&request);
}



void
InitFinalize::unregister_request(MPI_Request &request)
{
  Assert(requests.find(&request) != requests.end(),
         ExcMessage(
           "You tried to call unregister_request() with an invalid request."));

  requests.erase(&request);
}



std::set<MPI_Request *> InitFinalize::requests;



void
InitFinalize::finalize()
{
  if (!is_finalized)
    {
      // First, call the at_mpi_finalize() signal handler.
      signals.at_mpi_finalize();

      // make memory pool release all PETSc/Trilinos/MPI-based vectors that
      // are no longer used at this point. this is relevant because the static
      // object destructors run for these vectors at the end of the program
      // would run after MPI_Finalize is called, leading to errors

#ifdef DEAL_II_WITH_MPI
      // Before exiting, wait for nonblocking communication to complete:
      for (auto *request : requests)
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
      if (static_cast<bool>(libraries & InitializeLibrary::SLEPc) &&
          (finalize_petscslepc))
        {
          PetscErrorCode ierr = SlepcFinalize();
          AssertThrow(ierr == 0,
                      SLEPcWrappers::SolverBase::ExcSLEPcError(ierr));
        }
#  else
      // or just end PETSc if we did so
      if (static_cast<bool>(libraries & InitializeLibrary::PETSc) &&
          (finalize_petscslepc))
        {
          PetscErrorCode ierr = PetscFinalize();
          AssertThrow(ierr == 0, ExcPETScError(ierr));
        }
#  endif
#endif

#ifdef DEAL_II_WITH_P4EST
      // now end p4est and libsc
      // Note: p4est has no finalize function
      if (static_cast<bool>(libraries & InitializeLibrary::P4EST))
        sc_finalize();
#endif


      // Finalize Kokkos
      if (static_cast<bool>(libraries & InitializeLibrary::Kokkos))
        Kokkos::finalize();

        // only MPI_Finalize if we are running with MPI. We also need to do this
        // when running PETSc, because we initialize MPI ourselves before
        // calling PetscInitialize
#ifdef DEAL_II_WITH_MPI
      int       MPI_has_been_started = 0;
      const int ierr                 = MPI_Initialized(&MPI_has_been_started);
      AssertThrowMPI(ierr);
      if (static_cast<bool>(libraries & InitializeLibrary::MPI) &&
          (MPI_has_been_started))
        {
          if (std::uncaught_exceptions() > 0)
            {
              // do not try to call MPI_Finalize to avoid a deadlock.
            }
          else
            {
              const int ierr = MPI_Finalize();
              AssertNothrow(ierr == MPI_SUCCESS, dealii::ExcMPI(ierr));
            }
        }
#endif
      is_finalized = true;
    }
}



InitFinalize::~InitFinalize()
{
  finalize();
}


DEAL_II_NAMESPACE_CLOSE
