// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/base/multithread_info.h>

#include <cstddef>
#include <iostream>

#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#    include <deal.II/lac/vector_memory.h>
#    include <deal.II/lac/trilinos_vector.h>
#    include <deal.II/lac/trilinos_block_vector.h>
#  endif
#endif

#ifdef DEAL_II_WITH_PETSC
#  ifdef DEAL_II_WITH_MPI
#    include <petscsys.h>
#    include <deal.II/lac/petsc_block_vector.h>
#    include <deal.II/lac/petsc_parallel_block_vector.h>
#    include <deal.II/lac/petsc_vector.h>
#    include <deal.II/lac/petsc_parallel_vector.h>
#  endif
#endif

#ifdef DEAL_II_WITH_SLEPC
#  ifdef DEAL_II_WITH_MPI
#    include <slepcsys.h>
#  endif
#endif

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{

  namespace MPI
  {
#ifdef DEAL_II_WITH_MPI
    // Unfortunately, we have to work
    // around an oddity in the way PETSc
    // and some gcc versions interact. If
    // we use PETSc's MPI dummy
    // implementation, it expands the
    // calls to the two MPI functions
    // basically as ``(n_jobs=1, 0)'',
    // i.e. it assigns the number one to
    // the variable holding the number of
    // jobs, and then uses the comma
    // operator to let the entire
    // expression have the value zero. The
    // latter is important, since
    // ``MPI_Comm_size'' returns an error
    // code that we may want to check (we
    // don't here, but one could in
    // principle), and the trick with the
    // comma operator makes sure that both
    // the number of jobs is correctly
    // assigned, and the return value is
    // zero. Unfortunately, if some recent
    // versions of gcc detect that the
    // comma expression just stands by
    // itself, i.e. the result is not
    // assigned to another variable, then
    // they warn ``right-hand operand of
    // comma has no effect''. This
    // unwanted side effect can be
    // suppressed by casting the result of
    // the entire expression to type
    // ``void'' -- not beautiful, but
    // helps calming down unwarranted
    // compiler warnings...
    unsigned int n_mpi_processes (const MPI_Comm &mpi_communicator)
    {
      int n_jobs=1;
      (void) MPI_Comm_size (mpi_communicator, &n_jobs);

      return n_jobs;
    }


    unsigned int this_mpi_process (const MPI_Comm &mpi_communicator)
    {
      int rank=0;
      (void) MPI_Comm_rank (mpi_communicator, &rank);

      return rank;
    }


    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator)
    {
      MPI_Comm new_communicator;
      MPI_Comm_dup (mpi_communicator, &new_communicator);
      return new_communicator;
    }


    std::vector<unsigned int>
    compute_point_to_point_communication_pattern (const MPI_Comm &mpi_comm,
                                                  const std::vector<unsigned int> &destinations)
    {
      unsigned int myid = Utilities::MPI::this_mpi_process(mpi_comm);
      unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_comm);

      for (unsigned int i=0; i<destinations.size(); ++i)
        {
          Assert (destinations[i] < n_procs,
                  ExcIndexRange (destinations[i], 0, n_procs));
          Assert (destinations[i] != myid,
                  ExcMessage ("There is no point in communicating with ourselves."));
        }


      // let all processors
      // communicate the maximal
      // number of destinations they
      // have
      const unsigned int max_n_destinations
        = Utilities::MPI::max (destinations.size(), mpi_comm);

      // now that we know the number
      // of data packets every
      // processor wants to send, set
      // up a buffer with the maximal
      // size and copy our
      // destinations in there,
      // padded with -1's
      std::vector<unsigned int> my_destinations(max_n_destinations,
                                                numbers::invalid_unsigned_int);
      std::copy (destinations.begin(), destinations.end(),
                 my_destinations.begin());

      // now exchange these (we could
      // communicate less data if we
      // used MPI_Allgatherv, but
      // we'd have to communicate
      // my_n_destinations to all
      // processors in this case,
      // which is more expensive than
      // the reduction operation
      // above in MPI_Allreduce)
      std::vector<unsigned int> all_destinations (max_n_destinations * n_procs);
      MPI_Allgather (&my_destinations[0], max_n_destinations, MPI_UNSIGNED,
                     &all_destinations[0], max_n_destinations, MPI_UNSIGNED,
                     mpi_comm);

      // now we know who is going to
      // communicate with
      // whom. collect who is going
      // to communicate with us!
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
      // custom MIP_Op for
      // calculate_collective_mpi_min_max_avg
      void max_reduce ( const void *in_lhs_,
                        void *inout_rhs_,
                        int *len,
                        MPI_Datatype *)
      {
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
            // choose lower cpu index when tied to make operator cumutative
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
            // choose lower cpu index when tied to make operator cumutative
            if (inout_rhs->max_index > in_lhs->max_index)
              inout_rhs->max_index = in_lhs->max_index;
          }
      }
    }



    MinMaxAvg
    min_max_avg(const double my_value,
                const MPI_Comm &mpi_communicator)
    {
      MinMaxAvg result;

      const unsigned int my_id
        = dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
      const unsigned int numproc
        = dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);

      MPI_Op op;
      int ierr = MPI_Op_create((MPI_User_function *)&max_reduce, true, &op);
      AssertThrow(ierr == MPI_SUCCESS, ExcInternalError());

      MinMaxAvg in;
      in.sum = in.min = in.max = my_value;
      in.min_index = in.max_index = my_id;

      MPI_Datatype type;
      int lengths[]= {3,2};
      MPI_Aint displacements[]= {0,offsetof(MinMaxAvg, min_index)};
      MPI_Datatype types[]= {MPI_DOUBLE, MPI_INT};

      ierr = MPI_Type_struct(2, lengths, displacements, types, &type);
      AssertThrow(ierr == MPI_SUCCESS, ExcInternalError());

      ierr = MPI_Type_commit(&type);
      ierr = MPI_Allreduce (&in, &result, 1, type, op, mpi_communicator);
      AssertThrow(ierr == MPI_SUCCESS, ExcInternalError());

      ierr = MPI_Type_free (&type);
      AssertThrow(ierr == MPI_SUCCESS, ExcInternalError());

      ierr = MPI_Op_free(&op);
      AssertThrow(ierr == MPI_SUCCESS, ExcInternalError());

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
      :
      owns_mpi (true)
    {
      do_init(argc, argv);

      if (max_num_threads != numbers::invalid_unsigned_int)
        {
          // set maximum number of threads (also respecting the environment
          // variable that the called function evaluates) based on what
          // the user asked
          multithread_info.set_thread_limit(max_num_threads);
        }
      else
        // user wants automatic choice
        {
#ifdef DEAL_II_WITH_MPI
          // we need to figure out how many MPI processes there
          // are on the current node, as well as how many CPU cores
          // we have. for the first task, check what get_hostname()
          // returns and then to an allgather so each processor
          // gets the answer
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
          MPI_Allgather (&hostname_array[0], max_hostname_size, MPI_CHAR,
                         &all_hostnames[0], max_hostname_size, MPI_CHAR,
                         MPI_COMM_WORLD);

          // search how often our own hostname appears and the
          // how-manyth instance the current process represents
          unsigned int n_local_processes=0;
          unsigned int nth_process_on_host = 0;
          for (unsigned int i=0; i<MPI::n_mpi_processes(MPI_COMM_WORLD); ++i)
            if (std::string (&all_hostnames[0] + i*max_hostname_size) == hostname)
              {
                ++n_local_processes;
                if (i <= MPI::this_mpi_process (MPI_COMM_WORLD))
                  ++nth_process_on_host;
              }
          Assert (nth_process_on_host > 0, ExcInternalError());


          // compute how many cores each process gets. if the number does
          // not divide evenly, then we get one more core if we are
          // among the first few processes
          //
          // if the number would be zero, round up to one since every
          // process needs to have at least one thread
          const unsigned int n_threads
            = std::max(multithread_info.n_cpus / n_local_processes
                       +
                       (nth_process_on_host <= multithread_info.n_cpus % n_local_processes
                        ?
                        1
                        :
                        0),
                       1U);
#else
          const unsigned int n_threads = multithread_info.n_cpus;
#endif

          // finally set this number of threads
          multithread_info.set_thread_limit(n_threads);
        }
    }




    MPI_InitFinalize::MPI_InitFinalize (int    &argc,
                                        char ** &argv)
      :
      owns_mpi (true)
    {
      do_init(argc, argv);

      // set maximum number of threads (also respecting the environment
      // variable that the called function evaluates)
      multithread_info.set_thread_limit(1);
    }


    void
    MPI_InitFinalize::do_init(int    &argc,
                              char ** &argv)
    {
      static bool constructor_has_already_run = false;
      Assert (constructor_has_already_run == false,
              ExcMessage ("You can only create a single object of this class "
                          "in a program since it initializes the MPI system."));



#ifdef DEAL_II_WITH_MPI
      // if we have PETSc, we will initialize it and let it handle MPI.
      // Otherwise, we will do it.
      int MPI_has_been_started = 0;
      MPI_Initialized(&MPI_has_been_started);
      AssertThrow (MPI_has_been_started == 0,
                   ExcMessage ("MPI error. You can only start MPI once!"));

      int mpi_err, provided;
      // this works likempi_err = MPI_Init (&argc, &argv); but tells MPI that
      // we might use several threads but never call two MPI functions at the
      // same time. For an explanation see on why we do this see
      // http://www.open-mpi.org/community/lists/users/2010/03/12244.php
      int wanted = MPI_THREAD_SERIALIZED;
      mpi_err = MPI_Init_thread(&argc, &argv, wanted, &provided);
      AssertThrow (mpi_err == 0,
                   ExcMessage ("MPI could not be initialized."));

      // disable for now because at least some implementations always return MPI_THREAD_SINGLE.
      //Assert(max_num_threads==1 || provided != MPI_THREAD_SINGLE,
      //    ExcMessage("MPI reports that we are not allowed to use multiple threads."));
#else
      // make sure the compiler doesn't warn
      // about these variables
      (void)argc;
      (void)argv;
      (void)owns_mpi;
#endif

      // we are allowed to call MPI_Init ourselves and PETScInitialize will
      // detect this. This allows us to use MPI_Init_thread instead.
#ifdef DEAL_II_WITH_PETSC
#  ifdef DEAL_II_WITH_SLEPC
      // Initialize SLEPc (with PETSc):
      SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#  else
      // or just initialize PETSc alone:
      PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#  endif
#endif

      constructor_has_already_run = true;
    }


    MPI_InitFinalize::~MPI_InitFinalize()
    {
      // make memory pool release all PETSc/Trilinos/MPI-based vectors that are no
      // longer used at this point. this is relevant because the
      // static object destructors run for these vectors at the end of
      // the program would run after MPI_Finalize is called, leading
      // to errors

#ifdef DEAL_II_WITH_MPI
      // Start with the deal.II MPI vectors (need to do this before finalizing
      // PETSc because it finalizes MPI).  Delete vectors from the pools:
      GrowingVectorMemory<parallel::distributed::Vector<double> >
      ::release_unused_memory ();
      GrowingVectorMemory<parallel::distributed::BlockVector<double> >
      ::release_unused_memory ();
      GrowingVectorMemory<parallel::distributed::Vector<float> >
      ::release_unused_memory ();
      GrowingVectorMemory<parallel::distributed::BlockVector<float> >
      ::release_unused_memory ();

      // Next with Trilinos:
#  if defined(DEAL_II_WITH_TRILINOS)
      GrowingVectorMemory<TrilinosWrappers::MPI::Vector>
      ::release_unused_memory ();
      GrowingVectorMemory<TrilinosWrappers::MPI::BlockVector>
      ::release_unused_memory ();
#  endif
#endif


      // Now deal with PETSc (with or without MPI). Only delete the vectors if finalize hasn't
      // been called yet, otherwise this will lead to errors.
#ifdef DEAL_II_WITH_PETSC
      if ((PetscInitializeCalled == PETSC_TRUE)
          &&
          (PetscFinalizeCalled == PETSC_FALSE))
        {
          GrowingVectorMemory<PETScWrappers::MPI::Vector>
          ::release_unused_memory ();
          GrowingVectorMemory<PETScWrappers::MPI::BlockVector>
          ::release_unused_memory ();
          GrowingVectorMemory<PETScWrappers::Vector>
          ::release_unused_memory ();
          GrowingVectorMemory<PETScWrappers::BlockVector>
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


      // only MPI_Finalize if we are running with MPI. We also need to do this
      // when running PETSc, because we initialize MPI ourselves before calling
      // PetscInitialize
#ifdef DEAL_II_WITH_MPI
      int mpi_err = 0;

      int MPI_has_been_started = 0;
      MPI_Initialized(&MPI_has_been_started);
      if (Utilities::System::job_supports_mpi() == true && owns_mpi == true &&
          MPI_has_been_started != 0)
        {
          if (std::uncaught_exception())
            {
              std::cerr << "ERROR: Uncaught exception in MPI_InitFinalize on proc "
                        << this_mpi_process(MPI_COMM_WORLD)
                        << ". Skipping MPI_Finalize() to avoid a deadlock."
                        << std::endl;
            }
          else
            mpi_err = MPI_Finalize();
        }


      AssertThrow (mpi_err == 0,
                   ExcMessage ("An error occurred while calling MPI_Finalize()"));
#endif
    }


  } // end of namespace MPI

} // end of namespace Utilities

DEAL_II_NAMESPACE_CLOSE
