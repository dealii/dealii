// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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


#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/base/multithread_info.h>

#include <iostream>

#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#    include <deal.II/lac/vector_memory.h>
#    include <deal.II/lac/trilinos_vector.h>
#    include <deal.II/lac/trilinos_parallel_block_vector.h>
#  endif
#endif

#ifdef DEAL_II_WITH_PETSC
#  include <petscsys.h>
#  include <deal.II/lac/petsc_parallel_block_vector.h>
#  include <deal.II/lac/petsc_parallel_vector.h>
#endif

#ifdef DEAL_II_WITH_SLEPC
#    include <slepcsys.h>
#    include <deal.II/lac/slepc_solver.h>
#endif

#ifdef DEAL_II_WITH_P4EST
#   include <p4est_bits.h>
#endif

#ifdef DEAL_II_TRILINOS_WITH_ZOLTAN
#   include <zoltan_cpp.h>
#endif

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{

  namespace MPI
  {
    void wait(const MPI_Comm &mpi_communicator)
    {
      // see https://www.open-mpi.org/faq/?category=debugging#serial-debuggers
      const unsigned int np = Utilities::MPI::n_mpi_processes(mpi_communicator);
      const unsigned int myid = Utilities::MPI::this_mpi_process(mpi_communicator);
      for (unsigned int p =0; p < np; ++p )
        {
          if (p==myid)
            std::cout << "PID " << ::getpid()
                      << " rank " << myid
                      << " on " << Utilities::System::get_hostname()
                      << " ready for attach." << std::endl;
#ifdef DEAL_II_WITH_MPI
          MPI_Barrier(mpi_communicator);
#endif
        }

      if (myid == 0)
        {
          char a;
          std::cout << "Type any character and press Enter/Return to start: " << std::flush;
          std::cin >> a;
        }

#ifdef DEAL_II_WITH_MPI
      MPI_Barrier(mpi_communicator);
#endif
    }

#ifdef DEAL_II_WITH_MPI
    unsigned int n_mpi_processes (const MPI_Comm &mpi_communicator)
    {
      int n_jobs=1;
      const int ierr = MPI_Comm_size (mpi_communicator, &n_jobs);
      AssertThrowMPI(ierr);

      return n_jobs;
    }


    unsigned int this_mpi_process (const MPI_Comm &mpi_communicator)
    {
      int rank=0;
      const int ierr = MPI_Comm_rank (mpi_communicator, &rank);
      AssertThrowMPI(ierr);

      return rank;
    }


    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator)
    {
      MPI_Comm new_communicator;
      const int ierr = MPI_Comm_dup (mpi_communicator, &new_communicator);
      AssertThrowMPI(ierr);
      return new_communicator;
    }


    std::vector<unsigned int>
    compute_point_to_point_communication_pattern (const MPI_Comm &mpi_comm,
                                                  const std::vector<unsigned int> &destinations)
    {
      const unsigned int myid = Utilities::MPI::this_mpi_process(mpi_comm);
      const unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_comm);

      for (unsigned int i=0; i<destinations.size(); ++i)
        {
          Assert (destinations[i] < n_procs,
                  ExcIndexRange (destinations[i], 0, n_procs));
          Assert (destinations[i] != myid,
                  ExcMessage ("There is no point in communicating with ourselves."));
        }


      // let all processors communicate the maximal number of destinations
      // they have
      const unsigned int max_n_destinations
        = Utilities::MPI::max (destinations.size(), mpi_comm);

      if (max_n_destinations==0)
        // all processes have nothing to send/receive:
        return std::vector<unsigned int>();

      // now that we know the number of data packets every processor wants to
      // send, set up a buffer with the maximal size and copy our destinations
      // in there, padded with -1's
      std::vector<unsigned int> my_destinations(max_n_destinations,
                                                numbers::invalid_unsigned_int);
      std::copy (destinations.begin(), destinations.end(),
                 my_destinations.begin());

      // now exchange these (we could communicate less data if we used
      // MPI_Allgatherv, but we'd have to communicate my_n_destinations to all
      // processors in this case, which is more expensive than the reduction
      // operation above in MPI_Allreduce)
      std::vector<unsigned int> all_destinations (max_n_destinations * n_procs);
      const int ierr = MPI_Allgather (my_destinations.data(), max_n_destinations, MPI_UNSIGNED,
                                      all_destinations.data(), max_n_destinations, MPI_UNSIGNED,
                                      mpi_comm);
      AssertThrowMPI(ierr);

      // now we know who is going to communicate with whom. collect who is
      // going to communicate with us!
      std::vector<unsigned int> origins;
      for (unsigned int i=0; i<n_procs; ++i)
        for (unsigned int j=0; j<max_n_destinations; ++j)
          if (all_destinations[i*max_n_destinations + j] == myid)
            origins.push_back (i);
          else if (all_destinations[i*max_n_destinations + j] ==
                   numbers::invalid_unsigned_int)
            break;

      return origins;
    }


    namespace
    {
      // custom MIP_Op for calculate_collective_mpi_min_max_avg
      void max_reduce ( const void *in_lhs_,
                        void *inout_rhs_,
                        int *len,
                        MPI_Datatype *)
      {
        (void)len;
        const MinMaxAvg *in_lhs = static_cast<const MinMaxAvg *>(in_lhs_);
        MinMaxAvg *inout_rhs = static_cast<MinMaxAvg *>(inout_rhs_);

        Assert(*len==1, ExcInternalError());

        inout_rhs->sum += in_lhs->sum;
        if (inout_rhs->min>in_lhs->min)
          {
            inout_rhs->min = in_lhs->min;
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
            inout_rhs->max = in_lhs->max;
            inout_rhs->max_index = in_lhs->max_index;
          }
        else if (inout_rhs->max == in_lhs->max)
          {
            // choose lower cpu index when tied to make operator commutative
            if (inout_rhs->max_index > in_lhs->max_index)
              inout_rhs->max_index = in_lhs->max_index;
          }
      }
    }



    MinMaxAvg
    min_max_avg(const double my_value,
                const MPI_Comm &mpi_communicator)
    {
      // If MPI was not started, we have a serial computation and cannot run
      // the other MPI commands
      if (job_supports_mpi() == false)
        {
          MinMaxAvg result;
          result.sum = my_value;
          result.avg = my_value;
          result.min = my_value;
          result.max = my_value;
          result.min_index = 0;
          result.max_index = 0;

          return result;
        }

      // To avoid uninitialized values on some MPI implementations, provide
      // result with a default value already...
      MinMaxAvg result = { 0., std::numeric_limits<double>::max(),
                           -std::numeric_limits<double>::max(), 0, 0, 0.
                         };

      const unsigned int my_id
        = dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
      const unsigned int numproc
        = dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);

      MPI_Op op;
      int ierr = MPI_Op_create((MPI_User_function *)&max_reduce, true, &op);
      AssertThrowMPI(ierr);

      MinMaxAvg in;
      in.sum = in.min = in.max = my_value;
      in.min_index = in.max_index = my_id;

      MPI_Datatype type;
      int lengths[]= {3,2};
      MPI_Aint displacements[]= {0,offsetof(MinMaxAvg, min_index)};
      MPI_Datatype types[]= {MPI_DOUBLE, MPI_INT};

      ierr = MPI_Type_struct(2, lengths, displacements, types, &type);
      AssertThrowMPI(ierr);

      ierr = MPI_Type_commit(&type);
      AssertThrowMPI(ierr);
      ierr = MPI_Allreduce (&in, &result, 1, type, op, mpi_communicator);
      AssertThrowMPI(ierr);

      ierr = MPI_Type_free (&type);
      AssertThrowMPI(ierr);

      ierr = MPI_Op_free(&op);
      AssertThrowMPI(ierr);

      result.avg = result.sum / numproc;

      return result;
    }

#else

    unsigned int n_mpi_processes (const MPI_Comm &)
    {
      return 1;
    }



    unsigned int this_mpi_process (const MPI_Comm &)
    {
      return 0;
    }


    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator)
    {
      return mpi_communicator;
    }



    MinMaxAvg
    min_max_avg(const double my_value,
                const MPI_Comm &)
    {
      MinMaxAvg result;

      result.sum = my_value;
      result.avg = my_value;
      result.min = my_value;
      result.max = my_value;
      result.min_index = 0;
      result.max_index = 0;

      return result;
    }

#endif



    MPI_InitFinalize::MPI_InitFinalize (int    &argc,
                                        char ** &argv,
                                        const unsigned int max_num_threads)
    {
      static bool constructor_has_already_run = false;
      (void)constructor_has_already_run;
      Assert (constructor_has_already_run == false,
              ExcMessage ("You can only create a single object of this class "
                          "in a program since it initializes the MPI system."));


      int ierr;
#ifdef DEAL_II_WITH_MPI
      // if we have PETSc, we will initialize it and let it handle MPI.
      // Otherwise, we will do it.
      int MPI_has_been_started = 0;
      ierr = MPI_Initialized(&MPI_has_been_started);
      AssertThrowMPI(ierr);
      AssertThrow (MPI_has_been_started == 0,
                   ExcMessage ("MPI error. You can only start MPI once!"));

      int provided;
      // this works like ierr = MPI_Init (&argc, &argv); but tells MPI that
      // we might use several threads but never call two MPI functions at the
      // same time. For an explanation see on why we do this see
      // http://www.open-mpi.org/community/lists/users/2010/03/12244.php
      int wanted = MPI_THREAD_SERIALIZED;
      ierr = MPI_Init_thread(&argc, &argv, wanted, &provided);
      AssertThrowMPI(ierr);

      // disable for now because at least some implementations always return
      // MPI_THREAD_SINGLE.
      //Assert(max_num_threads==1 || provided != MPI_THREAD_SINGLE,
      //    ExcMessage("MPI reports that we are not allowed to use multiple threads."));
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
      AssertThrow (ierr == 0, SLEPcWrappers::SolverBase::ExcSLEPcError(ierr));
#  else
      // or just initialize PETSc alone:
      ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
#  endif

      // Disable PETSc exception handling. This just prints a large wall
      // of text that is not particularly helpful for what we do:
      PetscPopSignalHandler();
#endif

      //Initialize zoltan
#ifdef DEAL_II_TRILINOS_WITH_ZOLTAN
      float version;
      Zoltan_Initialize(argc, argv, &version);
#endif

#ifdef DEAL_II_WITH_P4EST
      //Initialize p4est and libsc components
      sc_init(MPI_COMM_WORLD, 0, 0, nullptr, SC_LP_SILENT);
      p4est_init (nullptr, SC_LP_SILENT);
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
          // first task, check what get_hostname() returns and then to an
          // allgather so each processor gets the answer
          //
          // in calculating the length of the string, don't forget the
          // terminating \0 on C-style strings
          const std::string hostname = Utilities::System::get_hostname();
          const unsigned int max_hostname_size = Utilities::MPI::max (hostname.size()+1,
                                                                      MPI_COMM_WORLD);
          std::vector<char> hostname_array (max_hostname_size);
          std::copy (hostname.c_str(), hostname.c_str()+hostname.size()+1,
                     hostname_array.begin());

          std::vector<char> all_hostnames(max_hostname_size *
                                          MPI::n_mpi_processes(MPI_COMM_WORLD));
          const int ierr = MPI_Allgather (hostname_array.data(), max_hostname_size, MPI_CHAR,
                                          all_hostnames.data(), max_hostname_size, MPI_CHAR,
                                          MPI_COMM_WORLD);
          AssertThrowMPI(ierr);

          // search how often our own hostname appears and the how-manyth
          // instance the current process represents
          unsigned int n_local_processes=0;
          unsigned int nth_process_on_host = 0;
          for (unsigned int i=0; i<MPI::n_mpi_processes(MPI_COMM_WORLD); ++i)
            if (std::string (all_hostnames.data() + i*max_hostname_size) == hostname)
              {
                ++n_local_processes;
                if (i <= MPI::this_mpi_process (MPI_COMM_WORLD))
                  ++nth_process_on_host;
              }
          Assert (nth_process_on_host > 0, ExcInternalError());


          // compute how many cores each process gets. if the number does not
          // divide evenly, then we get one more core if we are among the
          // first few processes
          //
          // if the number would be zero, round up to one since every process
          // needs to have at least one thread
          const unsigned int n_threads
            = std::max(MultithreadInfo::n_cores() / n_local_processes
                       +
                       (nth_process_on_host <= MultithreadInfo::n_cores() % n_local_processes
                        ?
                        1
                        :
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
      GrowingVectorMemory<LinearAlgebra::distributed::Vector<double> >
      ::release_unused_memory ();
      GrowingVectorMemory<LinearAlgebra::distributed::BlockVector<double> >
      ::release_unused_memory ();
      GrowingVectorMemory<LinearAlgebra::distributed::Vector<float> >
      ::release_unused_memory ();
      GrowingVectorMemory<LinearAlgebra::distributed::BlockVector<float> >
      ::release_unused_memory ();

      // Next with Trilinos:
#  if defined(DEAL_II_WITH_TRILINOS)
      GrowingVectorMemory<TrilinosWrappers::MPI::Vector>
      ::release_unused_memory ();
      GrowingVectorMemory<TrilinosWrappers::MPI::BlockVector>
      ::release_unused_memory ();
#  endif
#endif


      // Now deal with PETSc (with or without MPI). Only delete the vectors if
      // finalize hasn't been called yet, otherwise this will lead to errors.
#ifdef DEAL_II_WITH_PETSC
      if ((PetscInitializeCalled == PETSC_TRUE)
          &&
          (PetscFinalizeCalled == PETSC_FALSE))
        {
          GrowingVectorMemory<PETScWrappers::MPI::Vector>
          ::release_unused_memory ();
          GrowingVectorMemory<PETScWrappers::MPI::BlockVector>
          ::release_unused_memory ();

#  ifdef DEAL_II_WITH_SLEPC
          // and now end SLEPc (with PETSc)
          SlepcFinalize();
#  else
          // or just end PETSc.
          PetscFinalize();
#  endif
        }
#endif

#ifdef DEAL_II_WITH_P4EST
      // now end p4est and libsc
      // Note: p4est has no finalize function
      sc_finalize ();
#endif


      // only MPI_Finalize if we are running with MPI. We also need to do this
      // when running PETSc, because we initialize MPI ourselves before
      // calling PetscInitialize
#ifdef DEAL_II_WITH_MPI
      if (job_supports_mpi() == true)
        {
          if (std::uncaught_exception())
            {
              std::cerr << "ERROR: Uncaught exception in MPI_InitFinalize on proc "
                        << this_mpi_process(MPI_COMM_WORLD)
                        << ". Skipping MPI_Finalize() to avoid a deadlock."
                        << std::endl;
            }
          else
            {
              const int ierr = MPI_Finalize();
              AssertThrowMPI(ierr);
            }
        }
#endif
    }



    bool job_supports_mpi ()
    {
#ifdef DEAL_II_WITH_MPI
      int MPI_has_been_started = 0;
      const int ierr = MPI_Initialized(&MPI_has_been_started);
      AssertThrowMPI(ierr);

      return (MPI_has_been_started > 0);
#else
      return false;
#endif
    }



#include "mpi.inst"
  } // end of namespace MPI
} // end of namespace Utilities

DEAL_II_NAMESPACE_CLOSE
