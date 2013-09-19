// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __tests_tests_h
#define __tests_tests_h

// common definitions used in all the tests

#include <deal.II/base/config.h>
#include <deal.II/base/job_identifier.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/thread_management.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>


// given the name of a file, copy it to deallog
// and then delete it
inline void cat_file(const char *filename)
{
  std::ifstream in(filename);
  Assert (in, dealii::ExcIO());

  while (in)
    {
      std::string s;
      std::getline(in, s);
      dealii::deallog.get_file_stream() << s << "\n";
    }
  in.close();

  std::remove (filename);
}


#ifdef DEAL_II_WITH_PETSC
#include <petscsys.h>

inline void check_petsc_allocations()
{
  PetscStageLog stageLog;
  PetscLogGetStageLog(&stageLog);

  // I don't quite understand petsc and it looks like stageLog->stageInfo->classLog->classInfo[i].id
  // is always -1, so we look it up in stageLog->classLog, make sure it has the same number of entries:
  Assert(stageLog->stageInfo->classLog->numClasses == stageLog->classLog->numClasses, dealii::ExcInternalError());

  bool errors = false;
  for (int i=0;i<stageLog->stageInfo->classLog->numClasses;++i)
    {
      if (stageLog->stageInfo->classLog->classInfo[i].destructions
	  != stageLog->stageInfo->classLog->classInfo[i].creations)
	{
	  errors = true;
	  std::cerr << "ERROR: PETSc objects leaking of type '"
	       <<  stageLog->classLog->classInfo[i].name << "'"
	       << " with " << stageLog->stageInfo->classLog->classInfo[i].creations
	       << " creations and only "
	       << stageLog->stageInfo->classLog->classInfo[i].destructions
	       << " destructions." << std::endl;
	  }
    }

  if (errors)
    throw dealii::ExcMessage("PETSc memory leak");
}
#endif


// implicitly use the deal.II namespace everywhere, without us having to say
// so in each and every testcase
using namespace dealii;

// Function for initialize deallog. Normally, it should be called at
// the beginning of main() like
//
// initlog(__FILE__);
//
// This will open the correct output file, divert log output there and
// switch off screen output. If screen output is desired, provide the
// optional second argument as 'true'.
std::string deallogname;
std::ofstream deallogfile;

inline
void
initlog(const char *filename, bool console=false)
{
  deallogname = JobIdentifier::base_name(filename) + std::string("/output");
  deallogfile.open(deallogname.c_str());
  deallog.attach(deallogfile);
  if (!console)
    deallog.depth_console(0);

//TODO: Remove this line and replace by test_mode()
  deallog.threshold_float(1.e-8);
}


// append the directory name with an output file name that indicates
// the number of MPI processes
inline std::string
output_file_for_mpi (const std::string &directory)
{
#ifdef DEAL_II_WITH_MPI
  return (directory + "/ncpu_" +
          Utilities::int_to_string (Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)) +
          "/output");
#else
  return (directory + "/ncpu_1/output");
#endif
}


inline
void
mpi_initlog(const char *filename, bool console=false)
{
#ifdef DEAL_II_WITH_MPI
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  if (myid == 0)
    {
      deallogname = output_file_for_mpi(JobIdentifier::base_name(filename));
      deallogfile.open(deallogname.c_str());
      deallog.attach(deallogfile);
      if (!console)
        deallog.depth_console(0);

//TODO: Remove this line and replace by test_mode()
      deallog.threshold_float(1.e-8);
    }
#else
  // can't use this function if not using MPI
  Assert (false, ExcInternalError());
#endif
}



/* helper class to include the deallogs of all processors
   on proc 0 */
class MPILogInitAll
{
public:
  MPILogInitAll(const char *filename, bool console=false)
    : m_filename(filename)
  {
#ifdef DEAL_II_WITH_MPI
    unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
    deallogname = output_file_for_mpi(JobIdentifier::base_name(filename));
    if (myid != 0)
      deallogname = deallogname + Utilities::int_to_string(myid);
    deallogfile.open(deallogname.c_str());
    deallog.attach(deallogfile);
    if (!console)
      deallog.depth_console(0);

//TODO: Remove this line and replace by test_mode()
    deallog.threshold_float(1.e-8);
    deallog.push(Utilities::int_to_string(myid));
#else
    // can't use this function if not using MPI
    Assert (false, ExcInternalError());
#endif
  }

  ~MPILogInitAll()
  {
#ifdef DEAL_II_WITH_MPI
    unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
    unsigned int nproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

    deallog.pop();

    if (myid!=0)
      {
        deallog.detach();
        deallogfile.close();
      }

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEAL_II_WITH_PETSC
    check_petsc_allocations();
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  
    if (myid==0)
      {
        for (unsigned int i=1; i<nproc; ++i)
          {
            std::string filename = output_file_for_mpi(JobIdentifier::base_name(m_filename.c_str()))
                                   + Utilities::int_to_string(i);
            std::ifstream in(filename.c_str());
            Assert (in, ExcIO());

            while (in)
              {
                std::string s;
                std::getline(in, s);
                deallog.get_file_stream() << s << "\n";
              }
            in.close();
            std::remove (filename.c_str());
          }
      }
#else
    // can't use this function if not using MPI
    Assert (false, ExcInternalError());
#endif
  }
private:

  std::string m_filename;

};






#ifndef DEAL_II_STACKTRACE_SWITCH
#define DEAL_II_STACKTRACE_SWITCH

// A structure and a variable that are used to make sure that we do
// not show a stacktrace in out testcases, since this would lead to
// trouble with automatic comparison of output against a fixed
// baseline
struct SwitchOffStacktrace
{
  SwitchOffStacktrace ()
  {
    deal_II_exceptions::suppress_stacktrace_in_exceptions ();
  }
} deal_II_stacktrace_dummy;

#endif

DEAL_II_NAMESPACE_OPEN
namespace internal
{

  namespace Vector
  {

    extern unsigned int minimum_parallel_grain_size;

  }

  namespace SparseMatrix
  {

    extern unsigned int minimum_parallel_grain_size;

  }

}

// A structure and a variable that are used to set grainsizes for parallel
// mode smaller than they would otherwise be. this is used to test that the
// parallel algorithms in lac/ work alright.
struct SetGrainSizes
{

  SetGrainSizes ()
  {

    internal::Vector::minimum_parallel_grain_size = 2;

    internal::SparseMatrix::minimum_parallel_grain_size = 2;

  }

}
set_grain_sizes;



// spawn a thread that terminates the program after a certain time
// given by the environment variable WALLTIME. this makes
// sure we don't let jobs that deadlock on some mutex hang around
// forever. note that this is orthogonal to using "ulimit" in
// Makefile.rules, which only affects CPU time and consequently works
// on infinite loops but not deadlocks
//
// the actual sleep time isn't all that important. the point is that
// we need to kill deadlocked threads at one point, whenever that is
#ifdef DEAL_II_WITH_THREADS

struct DeadlockKiller
{
private:
  static void nuke_it ()
  {
    char *env = std::getenv("WALLTIME");

    if (env != 0)
      {
        std::istringstream conv (env);
        int delay;
        conv >> delay;
        if (conv)
          {
            sleep (delay);
            std::cerr << "Time's up: Killing job after "
                      << delay
                      << " seconds because it overran its allowed walltime."
                      << std::endl;
            std::abort ();
          }
        else
          {
            std::cerr << "Invalid value for WALLTIME environment variable."
                      << std::endl;
            std::abort ();
          }
      }
    else
      // environment variable is not set, so assume infinite wait time
      // and simply quit this thread
      {
      }
  }

public:
  DeadlockKiller ()
  {
    dealii::Threads::new_thread (&nuke_it);
  }
};

#endif

DEAL_II_NAMESPACE_CLOSE




#endif // __tests_tests_h
