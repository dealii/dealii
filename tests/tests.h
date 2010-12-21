//----------------------------  tests.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2006, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tests.h  ---------------------------

#ifndef __tests_tests_h
#define __tests_tests_h

// common definitions used in all the tests

#include <base/config.h>
#include <base/job_identifier.h>
#include <base/logstream.h>
#include <base/exceptions.h>
#include <base/utilities.h>
#include <base/thread_management.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

// implicitly use the deal.II namespace everywhere, without us having to say
// so in each and every testcase
using namespace dealii;


// overload floating point output operators for LogStream so that small
// numbers below a certain threshold are simply printed as zeros. this removes
// a number of possible places where output may differ depending on platform,
// compiler options, etc, simply because round-off is different.
inline
LogStream & operator << (LogStream &logstream,
                         const float d)
{
  if (std::fabs (d) < 1e-8)
    logstream.
#ifdef DEAL_II_TEMPL_OP_DISAMBIGUATION_BUG
      template
#endif
      operator << <float> (0.);
  else
    logstream.
#ifdef DEAL_II_TEMPL_OP_DISAMBIGUATION_BUG
      template
#endif
      operator << <float> (d);
  return logstream;
}



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
#ifdef DEAL_II_USE_MT

struct DeadlockKiller
{
  private:
    static void nuke_it ()
      {
	char * env = std::getenv("WALLTIME");

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
					   // environment variable is
					   // not set, so assume
					   // infinite wait time and
					   // simply quit this thread
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


// append the directory name with an output file name that indicates
// the number of MPI processes
std::string output_file_for_mpi (const std::string &directory)
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  return (directory + "/ncpu_" +
	  Utilities::int_to_string (Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD)) +
	  "/output");
#else
  return (directory + "/ncpu_1/output");
#endif
}



#endif // __tests_tests_h
