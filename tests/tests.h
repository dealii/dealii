// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2022 by the deal.II authors
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

#ifndef dealii_tests_h
#define dealii_tests_h

// common definitions used in all the tests

#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/cuda.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/job_identifier.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
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
#  include <stdlib.h>

struct DisableWindowsDebugRuntimeDialog
{
  DisableWindowsDebugRuntimeDialog()
  {
    _set_abort_behavior(0, _WRITE_ABORT_MSG);
  }
} deal_II_windows_crt_dialog;
#endif

// Redefine Assert as AssertThrow to make sure that the code is tested similarly
// in Release mode and in Debug mode. clang-format makes sure that this file is
// included after all regular header files but before all the other local header
// files.
#undef Assert
#define Assert AssertThrow

// implicitly use the deal.II namespace everywhere, without us having to say
// so in each and every testcase
using namespace dealii;


// ------------------------- Utility functions used in tests ------------------

/**
 * Go through the input stream @p in and filter out binary data for the key @p key .
 * The filtered stream is returned in @p out.
 */
void
filter_out_xml_key(std::istream &in, const std::string &key, std::ostream &out)
{
  std::string       line;
  bool              found   = false;
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
              line  = line.substr(0, line.find(">", 0) + 1) + line.substr(pos);
              found = false;
            }
          else
            line = line.substr(0, line.find(">", 0) + 1);
          out << line << std::endl;
        }
      else if (line.find(closing) != std::string::npos)
        {
          found = false;
          // remove everything before "<"
          line = line.substr(line.find("<", 0));
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
#  include <deal.II/lac/petsc_vector_base.h>
PetscReal
get_real_assert_zero_imag(const PETScWrappers::internal::VectorReference &a)
{
  Assert(a.imag() == 0.0, ExcInternalError());
  return a.real();
}
#endif

template <typename number>
number
get_real_assert_zero_imag(const std::complex<number> &a)
{
  Assert(a.imag() == 0.0, ExcInternalError());
  return a.real();
}

template <typename number>
number
get_real_assert_zero_imag(const number &a)
{
  return a;
}


// Cygwin has a different implementation for rand() which causes many tests to
// fail. This here is a reimplementation that gives the same sequence of numbers
// as a program that uses rand() on a typical linux machine. we put this into a
// namespace to not conflict with stdlib
namespace Testing
{
  /**
   * This function defines how to deal with signed overflow, which is undefined
   * behavior otherwise, in sums. Since unsigned overflow is well-defined there
   * is no reason to resort to this function.
   * The way we want to define the overflow is the following:
   * The value after the maximal value is the minimal one and the value
   * before the minimal one is the maximal one. Hence, we have to distinguish
   * three cases:
   * 1. $a+b>max$: This can only happen if both @p a and @p b are positive. By
   *               adding $min-max-1$ we are mapping $max+n$ to $min+n-1$ for
   *               $n>1$.
   * 2. $a+b<min$: This can only happen if both @p a and @p b are negative. By
   *               adding $max-min+1$ we are mapping $min-n$ to $max-n+1$ for
   *               $n>1$.
   * 3. $min<=a+b<=max$: No overflow.
   */
  template <typename Number>
  Number
  nonoverflow_add(Number a, Number b)
  {
    constexpr Number max = std::numeric_limits<Number>::max();
    constexpr Number min = std::numeric_limits<Number>::min();
    if (b > 0 && a > max - b)
      return (min + a) + (b - max) - 1;
    if (b < 0 && a < min - b)
      return (max + a) + (b - min) + 1;
    return a + b;
  }

  int
  rand(const bool reseed = false, const int seed = 1)
  {
    static int  r[32];
    static int  k;
    static bool inited = false;

    if (!inited || reseed)
      {
        // srand treats a seed 0 as 1 for some reason
        r[0]          = (seed == 0) ? 1 : seed;
        long int word = r[0];

        for (int i = 1; i < 31; ++i)
          {
            // This does:
            //   r[i] = (16807 * r[i-1]) % 2147483647;
            // but avoids overflowing 31 bits.
            const long int hi = word / 127773;
            const long int lo = word % 127773;
            word              = 16807 * lo - 2836 * hi;
            if (word < 0)
              word += 2147483647;
            r[i] = word;
          }
        k = 31;
        for (int i = 31; i < 34; ++i)
          {
            r[k % 32] = r[(k + 32 - 31) % 32];
            k         = (k + 1) % 32;
          }

        for (int i = 34; i < 344; ++i)
          {
            r[k % 32] =
              nonoverflow_add(r[(k + 32 - 31) % 32], r[(k + 32 - 3) % 32]);
            k = (k + 1) % 32;
          }
        inited = true;
        if (reseed == true)
          return 0; // do not generate new no
      }

    r[k % 32] = nonoverflow_add(r[(k + 32 - 31) % 32], r[(k + 32 - 3) % 32]);
    int ret   = r[k % 32];
    k         = (k + 1) % 32;
    return static_cast<unsigned int>(ret) >> 1;
  }

  // reseed our random number generator
  void
  srand(const int seed)
  {
    rand(true, seed);
  }
} // namespace Testing



// Get a uniformly distributed random value between min and max
template <typename T = double>
T
random_value(const T &min = static_cast<T>(0), const T &max = static_cast<T>(1))
{
  return min + (max - min) *
                 (static_cast<T>(Testing::rand()) / static_cast<T>(RAND_MAX));
}



// Construct a uniformly distributed random point, with each coordinate
// between min and max
template <int dim>
inline Point<dim>
random_point(const double &min = 0.0, const double &max = 1.0)
{
  Assert(max >= min, ExcMessage("Make sure max>=min"));
  Point<dim> p;
  for (unsigned int i = 0; i < dim; ++i)
    p[i] = random_value(min, max);
  return p;
}



// Construct a uniformly distributed random box, with each coordinate
// between min and max
template <int dim>
inline BoundingBox<dim>
random_box(const double &min = 0.0, const double &max = 1.0)
{
  Assert(max >= min, ExcMessage("Make sure max>=min"));
  std::vector<Point<dim>> p = {random_point<dim>(min, max),
                               random_point<dim>(min, max)};
  return BoundingBox<dim>(p);
}



// given the name of a file, copy it to deallog
// and then delete it
void
cat_file(const char *filename)
{
  {
    std::ifstream in(filename);
    Assert(in, dealii::ExcIO());
    deallog.get_file_stream() << in.rdbuf() << "\n";
  }

  std::remove(filename);
}


/*
 * Some tests (notably base/thread*, base/task*) create output that
 * comes out in random order. To make the output of these tests comparable,
 * we need to sort them.
 *
 * This function does just that with the file given. All streams writing
 * to this should be closed when calling this function.
 */
void
sort_file_contents(const std::string &filename)
{
  int error = std::system(
    (std::string("LC_ALL=C sort ") + filename + " -o " + filename).c_str());
  AssertThrow(error == 0, ExcInternalError());
}


/*
 * simple ADLER32 checksum for a range of chars
 */
template <class IT>
unsigned int
checksum(const IT &begin, const IT &end)
{
  AssertThrow(sizeof(unsigned int) == 4, ExcInternalError());
  AssertThrow(sizeof(*begin) == 1, ExcInternalError());

  unsigned int a = 1;
  unsigned int b = 0;

  IT it = begin;

  while (it != end)
    {
      a = (a + static_cast<unsigned char>(*it)) % 65521;
      b = (a + b) % 65521;
      ++it;
    }

  return (b << 16) | a;
}



/*
 * Replace all occurrences of ' &' by '& ' from the given file to hopefully be
 * more compiler independent with respect to __PRETTY_FUNCTION__
 *
 * Also, while GCC prepends the name by "virtual " if the function is virtual,
 * Intel's ICC does not do that, so filter that out as well.
 */
std::string
unify_pretty_function(const std::string &text)
{
  std::string t = text;
  t             = Utilities::replace_in_string(t, " &", " & ");
  t             = Utilities::replace_in_string(t, " & ,", "&,");
  t             = Utilities::replace_in_string(t, " & )", "&)");
  t             = Utilities::replace_in_string(t, " & ", "& ");
  t             = Utilities::replace_in_string(t, "virtual ", "");
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

#define check_solver_within_range(SolverType_COMMAND,                \
                                  CONTROL_COMMAND,                   \
                                  MIN_ALLOWED,                       \
                                  MAX_ALLOWED)                       \
  {                                                                  \
    const unsigned int previous_depth = deallog.depth_file(0);       \
    try                                                              \
      {                                                              \
        SolverType_COMMAND;                                          \
      }                                                              \
    catch (SolverControl::NoConvergence & exc)                       \
      {}                                                             \
    deallog.depth_file(previous_depth);                              \
    const unsigned int steps = CONTROL_COMMAND;                      \
    if (steps >= MIN_ALLOWED && steps <= MAX_ALLOWED)                \
      {                                                              \
        deallog << "Solver stopped within " << MIN_ALLOWED << " - "  \
                << MAX_ALLOWED << " iterations" << std::endl;        \
      }                                                              \
    else                                                             \
      {                                                              \
        deallog << "Solver stopped after " << steps << " iterations" \
                << std::endl;                                        \
      }                                                              \
  }

/*
 * Allow a test program to define a number that is very small to a given
 * tolerance to be output as zero. This is used e.g. for the output of float
 * numbers where roundoff difference can make the error larger than what we
 * have set for numdiff (that is appropriate for double variables).
 */
template <typename Number>
Number
filter_out_small_numbers(const Number number, const double tolerance)
{
  if (std::abs(number) < tolerance)
    return Number();
  else
    return number;
}


// ---------------- Functions used in initializing subsystems -----------------


/*
 * If we run 64 tests at the same time on a 64-core system, and each of
 * them runs 64 threads, then we get astronomical loads. Limit concurrency
 * to a fixed (small) number of threads, independent of the core count. The
 * limit defaults to 3 and can be overridden by the environment variable
 * TEST_N_THREADS.
 */
inline unsigned int
testing_max_num_threads()
{
  const int default_n_threads = 3;

  if (const char *penv = std::getenv("TEST_N_THREADS"))
    try
      {
        const int n_threads = Utilities::string_to_int(std::string(penv));
        return n_threads > 0 ? n_threads : default_n_threads;
      }
    catch (...)
      {
        return default_n_threads;
      }
  else
    return default_n_threads;
}

struct LimitConcurrency
{
  LimitConcurrency()
  {
    MultithreadInfo::set_thread_limit(testing_max_num_threads());
  }
} limit_concurrency;



#ifdef DEAL_II_WITH_PETSC
#  include <petscsys.h>

namespace
{
  void
  check_petsc_allocations()
  {
#  if DEAL_II_PETSC_VERSION_GTE(3, 2, 0)
    PetscStageLog stageLog;
    PetscLogGetStageLog(&stageLog);

    // I don't quite understand petsc and it looks like
    // stageLog->stageInfo->classLog->classInfo[i].id is always -1, so we look
    // it up in stageLog->classLog, make sure it has the same number of entries:
    Assert(stageLog->stageInfo->classLog->numClasses ==
             stageLog->classLog->numClasses,
           dealii::ExcInternalError());

    bool errors = false;
    for (int i = 0; i < stageLog->stageInfo->classLog->numClasses; ++i)
      {
        if (stageLog->stageInfo->classLog->classInfo[i].destructions !=
            stageLog->stageInfo->classLog->classInfo[i].creations)
          {
            errors = true;
            std::cerr
              << "ERROR: PETSc objects leaking of type '"
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
#  endif
  }
} // namespace
#endif


// Function to initialize deallog. Normally, it should be called at
// the beginning of main() like
//
// initlog();
//
// This will open the correct output file, divert log output there and
// switch off screen output. If screen output is desired, provide the
// optional first argument as 'true'.
std::string   deallogname;
std::ofstream deallogfile;

void
initlog(const bool                    console = false,
        const std::ios_base::fmtflags flags   = std::ios::showpoint |
                                              std::ios::left)
{
  deallogname = "output";
  deallogfile.open(deallogname);
  deallog.attach(deallogfile, true, flags);
  deallog.depth_console(console ? 10 : 0);
}


inline void
mpi_initlog(const bool                    console = false,
            const std::ios_base::fmtflags flags   = std::ios::showpoint |
                                                  std::ios::left)
{
#ifdef DEAL_II_WITH_MPI
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (myid == 0)
    {
      deallogname = "output";
      deallogfile.open(deallogname.c_str());
      deallog.attach(deallogfile, true, flags);
      deallog.depth_console(console ? 10 : 0);
    }
#else
  (void)console;
  (void)flags;
  // can't use this function if not using MPI
  Assert(false, ExcInternalError());
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
  MPILogInitAll(const bool                    console = false,
                const std::ios_base::fmtflags flags   = std::ios::showpoint |
                                                      std::ios::left)
  {
#ifdef DEAL_II_WITH_MPI
    const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
#else
    constexpr unsigned int myid = 0;
#endif
    if (myid == 0)
      {
        if (!deallog.has_file())
          {
            deallogfile.open("output");
            deallog.attach(deallogfile, true, flags);
          }
      }
    else
      {
        deallogname = "output" + Utilities::int_to_string(myid);
        deallogfile.open(deallogname.c_str());
        deallog.attach(deallogfile, true, flags);
      }

    deallog.depth_console(console ? 10 : 0);

    deallog.push(Utilities::int_to_string(myid));
  }

  ~MPILogInitAll()
  {
    // pop the prefix for the MPI rank of the current process
    deallog.pop();

#ifdef DEAL_II_WITH_MPI
    const unsigned int myid  = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    const unsigned int nproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

    if (myid != 0)
      {
        deallog.detach();
        deallogfile.close();
      }

    MPI_Barrier(MPI_COMM_WORLD);

#  ifdef DEAL_II_WITH_PETSC
    check_petsc_allocations();
    MPI_Barrier(MPI_COMM_WORLD);
#  endif

    if (myid == 0)
      {
        for (unsigned int i = 1; i < nproc; ++i)
          {
            std::string filename = "output" + Utilities::int_to_string(i);
            cat_file(filename.c_str());
          }
      }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
};


#ifdef DEAL_II_WITH_CUDA
// By default, all the ranks will try to access the device 0.
// If we are running with MPI support it is better to address different graphic
// cards for different processes even if only one node is used. The choice below
// is based on the MPI process id.
// MPI needs to be initialized before using this function.
void
init_cuda(const bool use_mpi = false)
{
#  ifndef DEAL_II_WITH_MPI
  Assert(use_mpi == false, ExcInternalError());
#  endif
  const unsigned int my_id =
    use_mpi ? Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) : 0;
  int         n_devices       = 0;
  cudaError_t cuda_error_code = cudaGetDeviceCount(&n_devices);
  AssertCuda(cuda_error_code);
  const int device_id = my_id % n_devices;
  cuda_error_code     = cudaSetDevice(device_id);
  AssertCuda(cuda_error_code);

  // In principle, we should be able to distribute the load better by
  // choosing a random graphics card. For some reason, this produces timeouts
  // on the tester we use mainly for the CUDA tests so we don't use the
  // following optimization by default.

  /*
  # ifndef DEAL_II_WITH_MPI
    Assert(use_mpi == false, ExcInternalError());
  #  endif
    const unsigned int my_id =
      use_mpi ? Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) : 0;
    int         device_id       = 0;
    int         n_devices       = 0;
    cudaError_t cuda_error_code = cudaGetDeviceCount(&n_devices);
    AssertCuda(cuda_error_code);
    if (my_id == 0)
      {
        Testing::srand(std::time(nullptr));
        device_id = Testing::rand() % n_devices;
      }
  #  ifdef DEAL_II_WITH_MPI
    if (use_mpi)
      MPI_Bcast(&device_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #  endif
    device_id       = (device_id + my_id) % n_devices;
    cuda_error_code = cudaSetDevice(device_id);
    AssertCuda(cuda_error_code);
  */
}
#endif



/* Override the tbb assertion handler in order to print a stacktrace:*/

#ifdef TBB_DO_ASSERT

#  include <tbb/tbb_stddef.h>

DEAL_II_NAMESPACE_OPEN
namespace deal_II_exceptions
{
  namespace internals
  {
    extern bool show_stacktrace;
  }
} // namespace deal_II_exceptions
DEAL_II_NAMESPACE_CLOSE


void
new_tbb_assertion_handler(const char *file,
                          int         line,
                          const char *expr,
                          const char *comment)
{
  // Print out the original assertion message
  std::cerr << "TBB assertion:" << std::endl;
  std::cerr << "Assertion " << expr << " failed on line " << line << " of file "
            << file << std::endl;
  std::cerr << "Detailed description: " << comment << std::endl;

  // Reenable abort and stacktraces:
  deal_II_exceptions::internals::allow_abort_on_exception = true;
  deal_II_exceptions::internals::show_stacktrace          = true;

  // And abort with a deal.II exception:
  Assert(false, ExcMessage("TBB Exception, see above"));
}

struct SetTBBAssertionHandler
{
  SetTBBAssertionHandler()
  {
    ::tbb::set_assertion_handler(new_tbb_assertion_handler);
  }
} set_tbb_assertion_handler;

#endif /*TBB_DO_ASSERT*/


// ---------------------- Adjust global variables in deal.II ------------------


DEAL_II_NAMESPACE_OPEN
/*
 * Now, change some global behavior of deal.II and supporting libraries:
 */

/* Disable stack traces: */

struct SwitchOffStacktrace
{
  SwitchOffStacktrace()
  {
    deal_II_exceptions::suppress_stacktrace_in_exceptions();
  }
} deal_II_stacktrace_dummy;



/* Enable floating point exceptions in debug mode and if we have
   detected that they are usable: */

struct EnableFPE
{
  EnableFPE()
  {
#if defined(DEBUG) && defined(DEAL_II_HAVE_FP_EXCEPTIONS)
    // enable floating point exceptions
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
  }
} deal_II_enable_fpe;


/* Set grainsizes for parallel mode smaller than they would otherwise be.
 * This is used to test that the parallel algorithms in lac/ work alright:
 */

namespace internal
{
  namespace VectorImplementation
  {
    extern unsigned int minimum_parallel_grain_size;
  }
  namespace SparseMatrixImplementation
  {
    extern unsigned int minimum_parallel_grain_size;
  }
} // namespace internal

struct SetGrainSizes
{
  SetGrainSizes()
  {
    internal::VectorImplementation::minimum_parallel_grain_size       = 2;
    internal::SparseMatrixImplementation::minimum_parallel_grain_size = 2;
  }
} set_grain_sizes;

DEAL_II_NAMESPACE_CLOSE

template <class T>
LogStream &
operator<<(LogStream &out, const std::vector<T> &v)
{
  for (std::size_t i = 0; i < v.size(); ++i)
    out << v[i] << (i == v.size() - 1 ? "" : " ");
  return out;
}

#endif // dealii_tests_h
