// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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

#ifndef dealii_tests_h
#define dealii_tests_h

// common definitions used in all the tests

#include <deal.II/base/config.h>
#include <deal.II/base/job_identifier.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/point.h>
#include <deal.II/base/patterns.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

#if defined(DEBUG) && defined(DEAL_II_HAVE_FP_EXCEPTIONS)
#  include <cfenv>
#endif


// silence extra diagnostics in the testsuite
#ifdef DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#endif


#ifdef DEAL_II_MSVC
// Under windows tests will hang and show a debugging dialog box from the
// debug CRT if an exception is encountered. Disable this:
#include <stdlib.h>

struct DisableWindowsDebugRuntimeDialog
{
  DisableWindowsDebugRuntimeDialog ()
  {
    _set_abort_behavior( 0, _WRITE_ABORT_MSG);
  }
} deal_II_windows_crt_dialog;
#endif

// implicitly use the deal.II namespace everywhere, without us having to say
// so in each and every testcase
using namespace dealii;


// ------------------------------ Utility functions used in tests -----------------------

/**
 * Go through the input stream @p in and filter out binary data for the key @p key .
 * The filtered stream is returned in @p out.
 */
void filter_out_xml_key(std::istream &in, const std::string &key, std::ostream &out)
{
  std::string line;
  bool found = false;
  const std::string opening = "<" + key;
  const std::string closing = "</" + key;
  while (std::getline(in, line))
    {
      if (line.find(opening) != std::string::npos &&
          line.find("binary") != std::string::npos)
        {
          found = true;
          // remove everything after ">" but keep things after "</"
          const auto pos = line.find(closing);
          if (pos != std::string::npos)
            {
              line = line.substr(0, line.find(">", 0)+1) + line.substr(pos);
              found = false;
            }
          else
            line = line.substr(0, line.find(">", 0)+1);
          out << line << std::endl;
        }
      else if (line.find(closing) != std::string::npos)
        {
          found = false;
          // remove everything before "<"
          line = line.substr(line.find("<",0));
          out << line << std::endl;
        }
      else if (!found)
        out << line << std::endl;
    }
}

/**
 * A function to return real part of the number and check that
 * its imaginary part is zero.
 */
#ifdef DEAL_II_WITH_PETSC
#include <deal.II/lac/petsc_vector_base.h>
PetscReal get_real_assert_zero_imag(const PETScWrappers::internal::VectorReference &a)
{
  Assert (a.imag() == 0.0, ExcInternalError());
  return a.real();
}
#endif

template <typename number>
number get_real_assert_zero_imag(const std::complex<number> &a)
{
  Assert (a.imag() == 0.0, ExcInternalError());
  return a.real();
}

template <typename number>
number get_real_assert_zero_imag(const number &a)
{
  return a;
}


// Cygwin has a different implementation for rand() which causes many tests to fail.
// This here is a reimplementation that gives the same sequence of numbers as a program
// that uses rand() on a typical linux machine.
// we put this into a namespace to not conflict with stdlib
namespace Testing
{
  int rand(const bool reseed=false,
           const int seed=1)
  {
    static int r[32];
    static int k;
    static bool inited=false;
    if (!inited || reseed)
      {
        //srand treats a seed 0 as 1 for some reason
        r[0]=(seed==0)?1:seed;

        for (int i=1; i<31; i++)
          {
            r[i] = (16807LL * r[i-1]) % 2147483647;
            if (r[i] < 0)
              r[i] += 2147483647;
          }
        k=31;
        for (int i=31; i<34; i++)
          {
            r[k%32] = r[(k+32-31)%32];
            k=(k+1)%32;
          }

        for (int i=34; i<344; i++)
          {
            r[k%32] = r[(k+32-31)%32] + r[(k+32-3)%32];
            k=(k+1)%32;
          }
        inited=true;
        if (reseed==true)
          return 0;// do not generate new no
      }

    r[k%32] = r[(k+32-31)%32] + r[(k+32-3)%32];
    int ret = r[k%32];
    k=(k+1)%32;
    return (unsigned int)ret >> 1;
  }

// reseed our random number generator
  void srand(const int seed)
  {
    rand(true, seed);
  }
}



// Get a uniformly distributed random value between min and max
template<typename T=double>
T random_value(const T &min=static_cast<T>(0),
               const T &max=static_cast<T>(1))
{
  return min+(max-min)*(static_cast<T>(Testing::rand())/static_cast<T>(RAND_MAX));
}



// Construct a uniformly distributed random point, with each coordinate
// between min and max
template<int dim>
inline Point<dim> random_point(const double &min=0.0,
                               const double &max=1.0)
{
  Assert(max >= min, ExcMessage("Make sure max>=min"));
  Point<dim> p;
  for (unsigned int i=0; i<dim; ++i)
    p[i] = random_value(min, max);
  return p;
}



// given the name of a file, copy it to deallog
// and then delete it
void cat_file(const char *filename)
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


/*
 * Some tests (notably base/thread*, base/task*) create output that
 * comes out in random order. To make the output of these tests comparable,
 * we need to sort them.
 *
 * This function does just that with the file given. All streams writing
 * to this should be closed when calling this function.
 */
void sort_file_contents (const std::string &filename)
{
  int error = std::system ((std::string ("LC_ALL=C sort ") + filename + " -o " + filename).c_str());
  AssertThrow (error == 0, ExcInternalError());
}


/*
 * simple ADLER32 checksum for a range of chars
 */
template <class IT>
unsigned int checksum(const IT &begin, const IT &end)
{
  AssertThrow(sizeof(unsigned int)==4, ExcInternalError());
  AssertThrow(sizeof(*begin)==1, ExcInternalError());

  unsigned int a = 1;
  unsigned int b = 0;

  IT it = begin;

  while (it != end)
    {
      a = (a + (unsigned char)*it) % 65521;
      b = (a + b) % 65521;
      ++it;
    }

  return (b << 16) | a;
}



/*
 * Replace all occurences of ' &' by '& ' from the given file to hopefully be
 * more compiler independent with respect to __PRETTY_FUNCTION__
 *
 * Also, while GCC prepends the name by "virtual " if the function is virtual,
 * Intel's ICC does not do that, so filter that out as well.
 */
std::string unify_pretty_function (const std::string &text)
{
  std::string t=text;
  t=Utilities::replace_in_string(t, " &", " & ");
  t=Utilities::replace_in_string(t, " & ,", "&,");
  t=Utilities::replace_in_string(t, " & )", "&)");
  t=Utilities::replace_in_string(t, " & ", "& ");
  t=Utilities::replace_in_string(t, "virtual ", "");
  return t;
}


/*
 * Test that a solver converged within a certain range of iteration steps.
 *
 * SolverType_COMMAND is the command to issue, CONTROL_COMMAND a function call
 * that returns the number of iterations (castable to unsigned int), and
 * MIN_ALLOWED, MAX_ALLOWED is the inclusive range of allowed iteration
 * steps.
 */

#define check_solver_within_range(SolverType_COMMAND, CONTROL_COMMAND, MIN_ALLOWED, MAX_ALLOWED) \
  {                                                                              \
    const unsigned int previous_depth = deallog.depth_file(0);                   \
    try                                                                          \
      {                                                                                \
        SolverType_COMMAND;                                                          \
      }                                                                                \
    catch  (SolverControl::NoConvergence &exc)                                                    \
      {}                                                                               \
    deallog.depth_file(previous_depth);                                          \
    const unsigned int steps = CONTROL_COMMAND;                                  \
    if (steps >= MIN_ALLOWED && steps <= MAX_ALLOWED)                            \
      {                                                                          \
        deallog << "Solver stopped within " << MIN_ALLOWED << " - "              \
                << MAX_ALLOWED << " iterations" << std::endl;                    \
      }                                                                          \
    else                                                                         \
      {                                                                          \
        deallog << "Solver stopped after " << steps << " iterations"             \
                << std::endl;                                                    \
      }                                                                          \
  }

/*
 * Allow a test program to define a number that is very small to a given
 * tolerance to be output as zero. This is used e.g. for the output of float
 * numbers where roundoff difference can make the error larger than what we
 * have set for numdiff (that is appropriate for double variables).
 */
template <typename Number>
Number filter_out_small_numbers (const Number number, const double tolerance)
{
  if (std::abs(number) < tolerance)
    return Number();
  else
    return number;
}

// ------------------------------ Functions used in initializing subsystems -------------------


/*
 * If we run 64 tests at the same time on a 64-core system, and
 * each of them runs 64 threads, then we get astronomical loads.
 * Limit concurrency to a fixed (small) number of threads, independent
 * of the core count.
 */
inline unsigned int testing_max_num_threads()
{
  return 5;
}

struct LimitConcurrency
{
  LimitConcurrency ()
  {
    MultithreadInfo::set_thread_limit (testing_max_num_threads());
  }
} limit_concurrency;



#ifdef DEAL_II_WITH_PETSC
#include <petscsys.h>

namespace
{
  void check_petsc_allocations()
  {
#if DEAL_II_PETSC_VERSION_GTE(3, 2, 0)
    PetscStageLog stageLog;
    PetscLogGetStageLog(&stageLog);

    // I don't quite understand petsc and it looks like
    // stageLog->stageInfo->classLog->classInfo[i].id is always -1, so we look
    // it up in stageLog->classLog, make sure it has the same number of entries:
    Assert(stageLog->stageInfo->classLog->numClasses == stageLog->classLog->numClasses,
           dealii::ExcInternalError());

    bool errors = false;
    for (int i=0; i<stageLog->stageInfo->classLog->numClasses; ++i)
      {
        if (stageLog->stageInfo->classLog->classInfo[i].destructions !=
            stageLog->stageInfo->classLog->classInfo[i].creations)
          {
            errors = true;
            std::cerr << "ERROR: PETSc objects leaking of type '"
                      << stageLog->classLog->classInfo[i].name << "'"
                      << " with "
                      << stageLog->stageInfo->classLog->classInfo[i].creations
                      << " creations and only "
                      << stageLog->stageInfo->classLog->classInfo[i].destructions
                      << " destructions." << std::endl;
          }
      }

    if (errors)
      throw dealii::ExcMessage("PETSc memory leak");
#endif
  }
}
#endif


// Function to initialize deallog. Normally, it should be called at
// the beginning of main() like
//
// initlog();
//
// This will open the correct output file, divert log output there and
// switch off screen output. If screen output is desired, provide the
// optional second argument as 'true'.
std::string deallogname;
std::ofstream deallogfile;

void
initlog(bool console=false)
{
  deallogname = "output";
  deallogfile.open(deallogname.c_str());
  deallog.attach(deallogfile);
  deallog.depth_console(console?10:0);
}


inline
void
mpi_initlog(bool console=false)
{
#ifdef DEAL_II_WITH_MPI
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  if (myid == 0)
    {
      deallogname = "output";
      deallogfile.open(deallogname.c_str());
      deallog.attach(deallogfile);
      deallog.depth_console(console?10:0);
    }
#else
  (void)console;
  // can't use this function if not using MPI
  Assert (false, ExcInternalError());
#endif
}



/**
 * A helper class that gives each MPI process its own output file
 * for the `deallog` stream, and at the end of the program (or,
 * more correctly, the end of the current object), concatenates them
 * all into the output file used on processor 0.
 */
struct MPILogInitAll
{
  MPILogInitAll(const bool console=false)
  {
#ifdef DEAL_II_WITH_MPI
    const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
    if (myid == 0)
      {
        if (!deallog.has_file())
          {
            deallogfile.open("output");
            deallog.attach(deallogfile);
          }
      }
    else
      {
        deallogname = "output" + Utilities::int_to_string(myid);
        deallogfile.open(deallogname.c_str());
        deallog.attach(deallogfile);
      }

    deallog.depth_console(console ? 10 : 0);

    deallog.push(Utilities::int_to_string(myid));
#else
    (void)console;
    // can't use this function if not using MPI
    Assert (false, ExcInternalError());
#endif
  }

  ~MPILogInitAll()
  {
#ifdef DEAL_II_WITH_MPI
    const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
    const unsigned int nproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

    // pop the prefix for the MPI rank of the current process
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
            std::string filename = "output" + Utilities::int_to_string(i);
            cat_file(filename.c_str());
          }
      }
    MPI_Barrier(MPI_COMM_WORLD);

#else
    // can't use this function if not using MPI
    Assert (false, ExcInternalError());
#endif
  }

};



/* Override the tbb assertion handler in order to print a stacktrace:*/

#ifdef TBB_DO_ASSERT

#include <tbb/tbb_stddef.h>

DEAL_II_NAMESPACE_OPEN
namespace deal_II_exceptions
{
  extern bool abort_on_exception;
  extern bool show_stacktrace;
}
DEAL_II_NAMESPACE_CLOSE


void new_tbb_assertion_handler(const char *file, int line, const char *expr,
                               const char *comment)
{
  // Print out the original assertion message
  std::cerr << "TBB assertion:" << std::endl;
  std::cerr << "Assertion " << expr << " failed on line " << line << " of file "
            << file << std::endl;
  std::cerr << "Detailed description: " << comment << std::endl;

  // Reenable abort and stacktraces:
  deal_II_exceptions::abort_on_exception = true;
  deal_II_exceptions::show_stacktrace = true;

  // And abort with a deal.II exception:
  Assert(false, ExcMessage("TBB Exception, see above"));
}

struct SetTBBAssertionHandler
{
  SetTBBAssertionHandler ()
  {
    ::tbb::set_assertion_handler(new_tbb_assertion_handler);
  }
} set_tbb_assertion_handler;

#endif /*TBB_DO_ASSERT*/


// ------------------------------ Adjust global variables in deal.II -----------------------


DEAL_II_NAMESPACE_OPEN
/*
 * Now, change some global behavior of deal.II and supporting libraries:
 */

/* Disable stack traces: */

struct SwitchOffStacktrace
{
  SwitchOffStacktrace ()
  {
    deal_II_exceptions::suppress_stacktrace_in_exceptions ();
  }
} deal_II_stacktrace_dummy;



/* Enable floating point exceptions in debug mode and if we have
   detected that they are usable: */

struct EnableFPE
{
  EnableFPE ()
  {
#if defined(DEBUG) && defined(DEAL_II_HAVE_FP_EXCEPTIONS)
    // enable floating point exceptions
    feenableexcept(FE_DIVBYZERO|FE_INVALID);
#endif
  }
} deal_II_enable_fpe;


/* Set grainsizes for parallel mode smaller than they would otherwise be.
 * This is used to test that the parallel algorithms in lac/ work alright:
 */

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

struct SetGrainSizes
{
  SetGrainSizes ()
  {
    internal::Vector::minimum_parallel_grain_size = 2;
    internal::SparseMatrix::minimum_parallel_grain_size = 2;
  }
} set_grain_sizes;

DEAL_II_NAMESPACE_CLOSE

/*
 * Do not use a template here to work around an overload resolution issue with clang and
 * enabled  C++11 mode.
 *
 * - Maier 2013
 */
LogStream &
operator << (LogStream &out,
             const std::vector<unsigned int> &v)
{
  for (unsigned int i=0; i<v.size(); ++i)
    out << v[i] << (i == v.size()-1 ? "" : " ");
  return out;
}

LogStream &
operator << (LogStream &out,
             const std::vector<long long unsigned int> &v)
{
  for (unsigned int i=0; i<v.size(); ++i)
    out << v[i] << (i == v.size()-1 ? "" : " ");
  return out;
}

LogStream &
operator << (LogStream &out,
             const std::vector<double> &v)
{
  for (unsigned int i=0; i<v.size(); ++i)
    out << v[i] << (i == v.size()-1 ? "" : " ");
  return out;
}


#endif // dealii_tests_h
