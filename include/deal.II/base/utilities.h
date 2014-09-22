// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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

#ifndef __deal2__utilities_h
#define __deal2__utilities_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>

#include <vector>
#include <utility>
#include <functional>
#include <string>

#ifdef DEAL_II_WITH_TRILINOS
#  include <Epetra_Comm.h>
#  include <Epetra_Map.h>
#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif
#endif

DEAL_II_NAMESPACE_OPEN


/**
 * A namespace for utility functions that are not particularly specific to
 * finite element computing or numerical programs, but nevertheless are needed
 * in various contexts when writing applications.
 *
 * @ingroup utilities
 * @author Wolfgang Bangerth, 2005
 */
namespace Utilities
{

  /**
   * Convert a number @p i to a string, with
   * as many digits as given to fill with
   * leading zeros.
   *
   * If the second parameter is left at its
   * default value, the number is not padded
   * with leading zeros. The result is then
   * the same as of the standard C function
   * <code>itoa()</code> had been called.
   */
  std::string
  int_to_string (const unsigned int i,
                 const unsigned int digits = numbers::invalid_unsigned_int);

  /**
   * Determine how many digits are needed to
   * represent numbers at most as large as
   * the given number.
   */
  unsigned int
  needed_digits (const unsigned int max_number);

  /**
   * Given a string, convert it to an
   * integer. Throw an assertion if that is
   * not possible.
   */
  int
  string_to_int (const std::string &s);

  /**
   * Return a string describing the dimensions of the object. Often,
   * functions in the deal.II library as well as in user codes need to
   * define a string containing the template dimensions of some
   * objects defined using two template parameters: dim (the
   * topological dimension of the object) and spacedim (the dimension
   * of the embedding Euclidean space).  Since in all deal.II classes,
   * by default spacedim is equal to dimension, the above string is
   * usually contracted to <dim>, instead of <dim,spacedim>. This
   * function returns a string containing "dim" if dim is equal to
   * spacedim, otherwhise it returns "dim,spacedim".
   */
  std::string dim_string(const int dim, const int spacedim);

  /**
   * Given a list of strings, convert it to a
   * list of integers. Throw an assertion if
   * that is not possible.
   */
  std::vector<int>
  string_to_int (const std::vector<std::string> &s);

  /**
   * Given a string, convert it to an
   * double. Throw an assertion if that is
   * not possible.
   */
  double
  string_to_double (const std::string &s);


  /**
   * Given a list of strings, convert it to a
   * list of doubles. Throw an assertion if
   * that is not possible.
   */
  std::vector<double>
  string_to_double (const std::vector<std::string> &s);

  /**
   * Given a string that contains text
   * separated by a @p delimiter, split it into
   * its components; for each component,
   * remove leading and trailing spaces.
   *
   * The default value of the delimiter is a
   * comma, so that the function splits comma
   * separated lists of strings.
   */
  std::vector<std::string>
  split_string_list (const std::string &s,
                     const char         delimiter = ',');

  /**
   * Take a text, usually a documentation or
   * something, and try to break it into
   * individual lines of text at most @p
   * width characters wide, by breaking at
   * positions marked by @p delimiter in the
   * text. If this is not possible, return
   * the shortest lines that are longer than
   * @p width.  The default value of the
   * delimiter is a space character. If
   * original_text contains newline
   * characters (\n), the string is split at
   * these locations, too.
   */
  std::vector<std::string>
  break_text_into_lines (const std::string &original_text,
                         const unsigned int width,
                         const char delimiter = ' ');

  /**
   * Return true if the given pattern
   * string appears in the first
   * position of the string.
   */
  bool
  match_at_string_start (const std::string &name,
                         const std::string &pattern);

  /**
   * Read a (signed) integer starting
   * at the position in @p name
   * indicated by the second
   * argument, and retun this integer
   * as a pair together with how many
   * characters it takes up in the
   * string.
   *
   * If no integer can be read at the
   * indicated position, return
   * (-1,numbers::invalid_unsigned_int)
   */
  std::pair<int, unsigned int>
  get_integer_at_position (const std::string &name,
                           const unsigned int position);

  /**
   * Generate a random number from a
   * normalized Gaussian probability
   * distribution centered around @p a and
   * with standard deviation @p sigma.
   *
   * This function is reentrant, i.e., it can safely be
   * called from multiple threads at the same time. However, if
   * so done, then there is no guarantee that each thread will
   * get the same sequence of numbers every time. Rather, the
   * produced sequence of random numbers will be apportioned to
   * the different threads in non-deterministic ways. If this
   * is a problem, for example for exactly reproducibility, then
   * you need to use separate random number facilities for separate
   * threads, rather than this global function. For example, the C++11
   * standard offers such objects, as does BOOST.
   *
   * @note Like the system function rand(), this function produces
   * the same sequence of random numbers every time a program is
   * started. This is an important property for debugging codes,
   * but it makes it impossible to really verify statistics
   * properties of a code. For rand(), you can call srand() to
   * "seed" the random number generator to get different sequences
   * of random numbers every time a program is called. However, this
   * function does not allow seeding the random number generator.
   * If you need this, as above, use one of the C++ or BOOST
   * facilities.
   */
  double
  generate_normal_random_number (const double a,
                                 const double sigma);


  /**
   * Calculate a fixed power, provided as a
   * template argument, of a number.
   *
   * This function provides an efficient way
   * to calculate things like
   * <code>t^N</code> where <code>N</code> is
   * a known number at compile time.
   *
   * Use this function as in
   * <code>fixed_power@<dim@> (n)</code>.
   */
  template <int N, typename T>
  T
  fixed_power (const T t);

  /**
   * Calculate a fixed power of an integer
   * number by a template expression where
   * both the number <code>a</code> and the
   * power <code>N</code> are compile-time
   * constants. This computes the result of
   * the power operation at compile time,
   * enabling its use e.g. in other
   * templates.
   *
   * Use this class as in
   * <code>fixed_int_power@<5,2@>::%value</code>
   * to compute 5<sup>2</sup>.
   */
  template <int a, int N>
  struct fixed_int_power
  {
    static const int value = a *fixed_int_power<a,N-1>::value;
  };

  /**
   * Base case for the power operation with
   * <code>N=0</code>, which gives the result
   * 1.
   */
  template <int a>
  struct fixed_int_power<a,0>
  {
    static const int value = 1;
  };

  /**
   * Optimized replacement for
   * <tt>std::lower_bound</tt> for
   * searching within the range of
   * column indices. Slashes
   * execution time by
   * approximately one half for the
   * present application, partly
   * because because the
   * binary search is replaced by a
   * linear search for small loop
   * lengths.
   *
   * Another reason for this function is
   * rather obscure:
   * when using the GCC libstdc++
   * function std::lower_bound,
   * complexity is O(log(N)) as
   * required. However, when using the
   * debug version of the GCC libstdc++
   * as we do when running the testsuite,
   * then std::lower_bound tests whether
   * the sequence is in fact partitioned
   * with respect to the pivot 'value'
   * (i.e. in essence that the sequence
   * is sorted as required for binary
   * search to work). However, verifying
   * this means that the complexity of
   * std::lower_bound jumps to O(N); we
   * call this function O(N) times below,
   * making the overall complexity
   * O(N**2). The consequence is that a
   * few tests with big meshes completely
   * run off the wall time limit for
   * tests and fail with the libstdc++
   * debug mode
   *
   * This function simply makes the
   * assumption that the sequence is
   * sorted, and we simply don't do the
   * additional check.
   */
  template<typename Iterator, typename T>
  Iterator
  lower_bound (Iterator  first,
               Iterator  last,
               const T  &val);


  /**
   * The same function as above, but taking
   * an argument that is used to compare
   * individual elements of the sequence of
   * objects pointed to by the iterators.
   */
  template<typename Iterator, typename T, typename Comp>
  Iterator
  lower_bound (Iterator   first,
               Iterator   last,
               const T   &val,
               const Comp comp);

  /**
   * Given a permutation vector (i.e. a
   * vector $p_0\ldots p_{N-1}$ where each
   * $p_i\in [0,N)$ and $p_i\neq p_j$ for
   * $i\neq j$), produce the reverse
   * permutation $q_i=N-1-p_i$.
   */
  std::vector<unsigned int>
  reverse_permutation (const std::vector<unsigned int> &permutation);

  /**
   * Given a permutation vector (i.e. a
   * vector $p_0\ldots p_{N-1}$ where each
   * $p_i\in [0,N)$ and $p_i\neq p_j$ for
   * $i\neq j$), produce the inverse
   * permutation $q_0\ldots q_{N-1}$ so that
   * $q_{p_i}=p_{q_i}=i$.
   */
  std::vector<unsigned int>
  invert_permutation (const std::vector<unsigned int> &permutation);

  /**
   * Given a permutation vector (i.e. a
   * vector $p_0\ldots p_{N-1}$ where each
   * $p_i\in [0,N)$ and $p_i\neq p_j$ for
   * $i\neq j$), produce the reverse
   * permutation $q_i=N-1-p_i$.
   */
  std::vector<unsigned long long int>
  reverse_permutation (const std::vector<unsigned long long int> &permutation);

  /**
   * Given a permutation vector (i.e. a
   * vector $p_0\ldots p_{N-1}$ where each
   * $p_i\in [0,N)$ and $p_i\neq p_j$ for
   * $i\neq j$), produce the inverse
   * permutation $q_0\ldots q_{N-1}$ so that
   * $q_{p_i}=p_{q_i}=i$.
   */
  std::vector<unsigned long long int>
  invert_permutation (const std::vector<unsigned long long int> &permutation);

  /**
   * A namespace for utility functions that
   * probe system properties.
   *
   * @ingroup utilities
   */
  namespace System
  {

    /**
     * Return the CPU load as returned by
     * "uptime". Note that the interpretation
     * of this number depends on the actual
     * number of processors in the
     * machine. This is presently only
     * implemented on Linux, using the
     * /proc/loadavg pseudo-file, on other
     * systems we simply return zero.
     */
    double get_cpu_load ();

    /**
     * Structure that hold information about
     * memory usage in kB. Used by
     * get_memory_stats(). See man 5 proc
     * entry /status for details.
     */
    struct MemoryStats
    {
      unsigned long int VmPeak; /** peak virtual memory size in kB */
      unsigned long int VmSize; /** current virtual memory size in kB */
      unsigned long int VmHWM; /** peak resident memory size in kB */
      unsigned long int VmRSS; /** current resident memory size in kB */
    };


    /**
     * Fills the @p stats structure with
     * information about the memory
     * consumption of this process. This is
     * only implemented on Linux.
     */
    void get_memory_stats (MemoryStats &stats);


    /**
     * Return the name of the host this
     * process runs on.
     */
    std::string get_hostname ();


    /**
     * Return the present time as HH:MM:SS.
     */
    std::string get_time ();

    /**
     * Return whether (i) deal.II has
     * been compiled to support MPI
     * (for example by compiling with
     * <code>CXX=mpiCC</code>) and if
     * so whether (ii)
     * <code>MPI_Init()</code> has
     * been called (for example using
     * the
     * Utilities::System::MPI_InitFinalize
     * class). In other words, the
     * result indicates whether the
     * current job is running under
     * MPI.
     *
     * @note The function does not
     * take into account whether an
     * MPI job actually runs on more
     * than one processor or is, in
     * fact, a single-node job that
     * happens to run under MPI.
     */
    bool job_supports_mpi ();

    /**
     * Alias for job_supports_mpi().
     *
     * @deprecated
     */
    bool program_uses_mpi () DEAL_II_DEPRECATED;

    /**
     * @name Functions that work
     * in parallel via MPI. The
     * functions following here
     * are all deprecated and have
     * been moved to namespace
     * Utilities::MPI.
     */
    /** @{ */

    /**
     * This function is an alias for
     * Utilities::MPI::n_mpi_processes.
     *
     * @deprecated
     */
    unsigned int get_n_mpi_processes (const MPI_Comm &mpi_communicator) DEAL_II_DEPRECATED;

    /**
     * This function is an alias for
     * Utilities::MPI::this_mpi_process.
     *
     * @deprecated
     */
    unsigned int get_this_mpi_process (const MPI_Comm &mpi_communicator) DEAL_II_DEPRECATED;


    /**
     * This function is an alias for
     * Utilities::MPI::compute_point_to_point_communication_pattern.
     *
     * @deprecated
     */
    using
    Utilities::MPI::compute_point_to_point_communication_pattern;


    /**
     * This function is an alias for
     * Utilities::MPI::duplicate_communicator.
     *
     * @deprecated
     */
    using Utilities::MPI::duplicate_communicator;

    /**
     * An alias for Utilities::MPI::MinMaxAvg.
     *
     * @deprecated
     */
    using Utilities::MPI::MinMaxAvg;

    /**
     * An alias for Utilities::MPI::min_max_avg.
     *
     * @deprecated
     */
    void
    calculate_collective_mpi_min_max_avg (const MPI_Comm &mpi_communicator,
                                          const double my_value,
                                          MinMaxAvg &result) DEAL_II_DEPRECATED;

    /**
     * An alias for Utilities::MPI::MPI_InitFinalize.
     *
     * @deprecated
     */
    using Utilities::MPI::MPI_InitFinalize;

    /** @} */
  }


#ifdef DEAL_II_WITH_TRILINOS
  /**
   * This namespace provides some of the basic structures used in the
   * initialization of the Trilinos objects (e.g., matrices, vectors, and
   * preconditioners).
   */
  namespace Trilinos
  {
    /**
     * Returns a Trilinos Epetra_Comm
     * object needed for creation of
     * Epetra_Maps.
     *
     * If deal.II has been configured to use
     * a compiler that does not support MPI
     * then the resulting communicator will
     * be a serial one. Otherwise, the
     * communicator will correspond to
     * MPI_COMM_WORLD, i.e. a communicator
     * that encompasses all processes within
     * this MPI universe.
     */
    const Epetra_Comm &comm_world();

    /**
     * Returns a Trilinos Epetra_Comm
     * object needed for creation of
     * Epetra_Maps.
     *
     * If deal.II has been configured to use
     * a compiler that does not support MPI
     * then the resulting communicator will
     * be a serial one. Otherwise, the
     * communicator will correspond to
     * MPI_COMM_SELF, i.e. a communicator
     * that comprises only this one
     * processor.
     */
    const Epetra_Comm &comm_self();

    /**
     * Given a communicator, duplicate it. If
     * the given communicator is serial, that
     * means to just return a copy of
     * itself. On the other hand, if it is
     * %parallel, we duplicate the underlying
     * MPI_Comm object: we create a separate
     * MPI communicator that contains the
     * same processors and in the same order
     * but has a separate identifier distinct
     * from the given communicator. The
     * function returns a pointer to a new
     * object of a class derived from
     * Epetra_Comm. The caller of this
     * function needs to assume ownership of
     * this function. The returned object
     * should be destroyed using the
     * destroy_communicator() function.
     *
     * This facility is used to separate
     * streams of communication. For example,
     * a program could simply use
     * MPI_Comm_World for everything. But it
     * is easy to come up with scenarios
     * where sometimes not all processors
     * participate in a communication that is
     * intended to be global -- for example
     * if we assemble a matrix on a coarse
     * mesh with fewer cells than there are
     * processors, some processors may not
     * sync their matrices with the rest
     * because they haven't written into it
     * because they own no cells. That's
     * clearly a bug. However, if these
     * processors just continue their work,
     * and the next %parallel operation
     * happens to be a sync on a different
     * matrix, then the sync could succeed --
     * by accident, since different
     * processors are talking about different
     * matrices.
     *
     * This kind of situation can be avoided
     * if we use different communicators for
     * different matrices which reduces the
     * likelihood that communications meant
     * to be separate aren't recognized as
     * such just because they happen on the
     * same communicator. In addition, it is
     * conceivable that some MPI operations
     * can be parallelized using multiple
     * threads because their communicators
     * identifies the communication in
     * question, not their relative timing as
     * is the case in a sequential program
     * that just uses a single communicator.
     */
    Epetra_Comm *
    duplicate_communicator (const Epetra_Comm &communicator);

    /**
     * Given an Epetra communicator that was
     * created by the
     * duplicate_communicator() function,
     * destroy the underlying MPI
     * communicator object and reset the
     * Epetra_Comm object to a the result of
     * comm_self().
     *
     * It is necessary to call this function
     * at the time when the result of
     * duplicate_communicator() is no longer
     * needed. The reason is that in that
     * function, we first create a new
     * MPI_Comm object and then create an
     * Epetra_Comm around it. While we can
     * take care of destroying the latter, it
     * doesn't destroy the communicator since
     * it can only assume that it may also be
     * still used by other objects in the
     * program. Consequently, we have to take
     * care of destroying it ourselves,
     * explicitly.
     *
     * This function does exactly
     * that. Because this has to happen while
     * the Epetra_Comm object is still
     * around, it first resets the latter and
     * then destroys the communicator object.
     *
     * @note If you call this function on an
     * Epetra_Comm object that is not created
     * by duplicate_communicator(), you are
     * likely doing something quite
     * wrong. Don't do this.
     */
    void
    destroy_communicator (Epetra_Comm &communicator);

    /**
     * Return the number of MPI processes
     * there exist in the given communicator
     * object. If this is a sequential job,
     * it returns 1.
     */
    unsigned int get_n_mpi_processes (const Epetra_Comm &mpi_communicator);

    /**
     * Return the number of the present MPI
     * process in the space of processes
     * described by the given
     * communicator. This will be a unique
     * value for each process between zero
     * and (less than) the number of all
     * processes (given by
     * get_n_mpi_processes()).
     */
    unsigned int get_this_mpi_process (const Epetra_Comm &mpi_communicator);

    /**
     * Given a Trilinos Epetra map, create a
     * new map that has the same subdivision
     * of elements to processors but uses the
     * given communicator object instead of
     * the one stored in the first
     * argument. In essence, this means that
     * we create a map that communicates
     * among the same processors in the same
     * way, but using a separate channel.
     *
     * This function is typically used with a
     * communicator that has been obtained by
     * the duplicate_communicator() function.
     */
    Epetra_Map
    duplicate_map (const Epetra_BlockMap  &map,
                   const Epetra_Comm &comm);
  }

#endif


}


// --------------------- inline functions

namespace Utilities
{
  template <int N, typename T>
  inline
  T fixed_power (const T n)
  {
    Assert (N>0, ExcNotImplemented());
    switch (N)
      {
      case 1:
        return n;
      case 2:
        return n*n;
      case 3:
        return n*n*n;
      case 4:
        return n*n*n*n;
      default:
        T result = n;
        for (int d=1; d<N; ++d)
          result *= n;
        return result;
      }
  }



  template<typename Iterator, typename T>
  inline
  Iterator
  lower_bound (Iterator  first,
               Iterator  last,
               const T  &val)
  {
    return Utilities::lower_bound (first, last, val,
                                   std::less<T>());
  }



  template<typename Iterator, typename T, typename Comp>
  inline
  Iterator
  lower_bound (Iterator    first,
               Iterator    last,
               const T    &val,
               const Comp  comp)
  {
    // verify that the two iterators are properly ordered. since
    // we need operator- for the iterator type anyway, do the
    // test as follows, rather than via 'last >= first'
    Assert (last - first >= 0,
            ExcMessage ("The given iterators do not satisfy the proper ordering."));

    unsigned int len = static_cast<unsigned int>(last-first);

    if (len==0)
      return first;

    while (true)
      {
        // if length equals 8 or less,
        // then do a rolled out
        // search. use a switch without
        // breaks for that and roll-out
        // the loop somehow
        if (len < 8)
          {
            switch (len)
              {
              case 7:
                if (!comp(*first, val))
                  return first;
                ++first;
              case 6:
                if (!comp(*first, val))
                  return first;
                ++first;
              case 5:
                if (!comp(*first, val))
                  return first;
                ++first;
              case 4:
                if (!comp(*first, val))
                  return first;
                ++first;
              case 3:
                if (!comp(*first, val))
                  return first;
                ++first;
              case 2:
                if (!comp(*first, val))
                  return first;
                ++first;
              case 1:
                if (!comp(*first, val))
                  return first;
                return first+1;
              default:
                // indices seem
                // to not be
                // sorted
                // correctly!? or
                // did len
                // become==0
                // somehow? that
                // shouln't have
                // happened
                Assert (false, ExcInternalError());
              }
          }



        const unsigned int half   = len >> 1;
        const Iterator     middle = first + half;

        // if the value is larger than
        // that pointed to by the
        // middle pointer, then the
        // insertion point must be
        // right of it
        if (comp(*middle, val))
          {
            first = middle + 1;
            len  -= half + 1;
          }
        else
          len = half;
      }
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
