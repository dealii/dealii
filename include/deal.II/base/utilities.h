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

#ifndef dealii_utilities_h
#define dealii_utilities_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <vector>
#include <utility>
#include <functional>
#include <string>
#include <tuple>
#include <type_traits>

#ifdef DEAL_II_WITH_TRILINOS
#  include <Epetra_Comm.h>
#  include <Epetra_Map.h>
#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif
#endif

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>

#ifdef DEAL_II_WITH_ZLIB
#  include <boost/iostreams/stream.hpp>
#  include <boost/iostreams/filtering_stream.hpp>
#  include <boost/iostreams/device/back_inserter.hpp>
#  include <boost/iostreams/filter/gzip.hpp>
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
   * Return a string of the form "deal.II version x.y.z" where "x.y.z"
   * identifies the version of deal.II you are using. This information
   * is also provided by the DEAL_II_PACKAGE_NAME and
   * DEAL_II_PACKAGE_VERSION preprocessor variables.
   */
  std::string dealii_version_string ();

  /**
   * Convert a number @p value to a string, with as many digits as given to
   * fill with leading zeros.
   *
   * If the second parameter is left at its default value, the number is not
   * padded with leading zeros. The result is then the same as if the standard
   * C function <code>itoa()</code> had been called.
   *
   * When calling this function signed integers are implicitly converted to
   * unsigned integers and long integers might experience an overflow.
   *
   * @note The use of this function is discouraged and users should use
   * <code>Utilities::to_string()</code> instead. In its current
   * implementation the function simply calls <code>to_string@<unsigned
   * int@>()</code>.
   */
  std::string
  int_to_string (const unsigned int value,
                 const unsigned int digits = numbers::invalid_unsigned_int);

  /**
   * Convert a number @p value to a string, with @p digits characters. The
   * string is padded with leading zeros, after a possible minus sign.
   * Therefore the total number of padding zeros is @p digits minus any signs,
   * decimal points and digits of @p value.
   *
   * If the second parameter is left at its default value, the number is not
   * padded with leading zeros. The result is then the same as if the boost
   * function <code>lexical_cast@<std::string@>()</code> had been called.
   */
  template <typename number>
  std::string
  to_string (const number value,
             const unsigned int digits = numbers::invalid_unsigned_int);

  /**
   * Determine how many digits are needed to represent numbers at most as
   * large as the given number.
   */
  unsigned int
  needed_digits (const unsigned int max_number);

  /**
   * Given a string, convert it to an integer. Throw an assertion if that is
   * not possible.
   */
  int
  string_to_int (const std::string &s);

  /**
   * Return a string describing the dimensions of the object. Often, functions
   * in the deal.II library as well as in user codes need to define a string
   * containing the template dimensions of some objects defined using two
   * template parameters: dim (the topological dimension of the object) and
   * spacedim (the dimension of the embedding Euclidean space).  Since in all
   * deal.II classes, by default spacedim is equal to dimension, the above
   * string is usually contracted to "<dim>", instead of "<dim,spacedim>".
   * This function returns a string containing "dim" if dim is equal to
   * spacedim, otherwise it returns "dim,spacedim".
   */
  std::string dim_string(const int dim, const int spacedim);

  /**
   * Given a list of strings, convert it to a list of integers. Throw an
   * assertion if that is not possible.
   */
  std::vector<int>
  string_to_int (const std::vector<std::string> &s);

  /**
   * Given a string, convert it to an double. Throw an assertion if that is
   * not possible.
   */
  double
  string_to_double (const std::string &s);


  /**
   * Given a list of strings, convert it to a list of doubles. Throw an
   * assertion if that is not possible.
   */
  std::vector<double>
  string_to_double (const std::vector<std::string> &s);


  /**
   * Given a string that contains text separated by a @p delimiter, split it
   * into its components; for each component, remove leading and trailing
   * spaces. The default value of the delimiter is a comma, so that the
   * function splits comma separated lists of strings.
   *
   * To make data input from tables simpler, if the input string ends in a
   * delimiter (possibly followed by an arbitrary amount of whitespace), then
   * this last delimiter is ignored. For example,
   * @code
   *   Utilities::split_string_list("abc; def; ghi; ", ';');
   * @endcode
   * yields the same 3-element list of output <code>{"abc","def","ghi"}</code>
   * as you would get if the input had been
   * @code
   *   Utilities::split_string_list("abc; def; ghi", ';');
   * @endcode
   * or
   * @code
   *   Utilities::split_string_list("abc; def; ghi;", ';');
   * @endcode
   * As a consequence of this rule, a call like
   * @code
   *   Utilities::split_string_list(" ; ", ';');
   * @endcode
   * yields a one-element list. Because of the trimming of whitespace, the
   * single element is the empty string.
   *
   * This function can digest the case that the delimiter is a space. In this
   * case, it returns all words in the string. Combined with the rules above,
   * this implies that
   * @code
   *   Utilities::split_string_list("abc def ghi ", ' ');
   * @endcode
   * yields again the 3-element list of output
   * <code>{"abc","def","ghi"}</code> from above despite the presence of space
   * at the end of the string. Furthermore,
   * @code
   *   Utilities::split_string_list("      ", ' ');
   * @endcode
   * yields an empty list regardless of the number of spaces in the string.
   */
  std::vector<std::string>
  split_string_list (const std::string &s,
                     const std::string &delimiter = ",");


  /**
   * Specialization of split_string_list() for the case where the delimiter
   * is a single char.
   */
  std::vector<std::string>
  split_string_list (const std::string &s,
                     const char delimiter);


  /**
   * Take a text, usually a documentation or something, and try to break it
   * into individual lines of text at most @p width characters wide, by
   * breaking at positions marked by @p delimiter in the text. If this is not
   * possible, return the shortest lines that are longer than @p width.  The
   * default value of the delimiter is a space character. If original_text
   * contains newline characters (\n), the string is split at these locations,
   * too.
   */
  std::vector<std::string>
  break_text_into_lines (const std::string &original_text,
                         const unsigned int width,
                         const char delimiter = ' ');

  /**
   * Return true if the given pattern string appears in the first position of
   * the string.
   */
  bool
  match_at_string_start (const std::string &name,
                         const std::string &pattern);

  /**
   * Read a (signed) integer starting at the position in @p name indicated by
   * the second argument, and return this integer as a pair together with how
   * many characters it takes up in the string.
   *
   * If no integer can be read at the indicated position, return
   * (-1,numbers::invalid_unsigned_int)
   */
  std::pair<int, unsigned int>
  get_integer_at_position (const std::string &name,
                           const unsigned int position);

  /**
   * Return a string with all occurrences of @p from in @p input replaced by
   * @p to.
   */
  std::string replace_in_string(const std::string &input,
                                const std::string &from,
                                const std::string &to);

  /**
   * Return a string with all standard whitespace characters (including
   * '<tt>\\t</tt>', '<tt>\\n</tt>', and '<tt>\\r</tt>') at the beginning and
   * end of @p input removed.
   */
  std::string
  trim(const std::string &input);

  /**
   * Generate a random number from a normalized Gaussian probability
   * distribution centered around @p a and with standard deviation @p sigma.
   *
   * This function is reentrant, i.e., it can safely be called from multiple
   * threads at the same time. In addition, each thread will get the same
   * sequence of numbers every time. On the other hand, if you run
   * Threads::Task objects via the Threading Building Blocks, then tasks will
   * be assigned to mostly random threads, and may get a different sequence of
   * random numbers in different runs of the program, since a previous task
   * may already have consumed the first few random numbers generated for the
   * thread you're on. If this is a problem, you need to create your own
   * random number generator objects every time you want to start from a
   * defined point.
   *
   * @note Like the system function rand(), this function produces the same
   * sequence of random numbers every time a program is started. This is an
   * important property for debugging codes, but it makes it impossible to
   * really verify statistics properties of a code. For rand(), you can call
   * srand() to "seed" the random number generator to get different sequences
   * of random numbers every time a program is called. However, this function
   * does not allow seeding the random number generator. If you need this, as
   * above, use one of the C++ or BOOST facilities.
   */
  double
  generate_normal_random_number (const double a,
                                 const double sigma);


  /**
   * Calculate a fixed power, provided as a template argument, of a number.
   *
   * This function provides an efficient way to calculate things like
   * <code>t^N</code> where <code>N</code> is a known number at compile time.
   *
   * Use this function as in <code>fixed_power@<dim@> (n)</code>.
   */
  template <int N, typename T>
  T
  fixed_power (const T t);

  /**
   * Calculate a fixed power of an integer number by a template expression
   * where both the number <code>a</code> and the power <code>N</code> are
   * compile-time constants. This computes the result of the power operation
   * at compile time, enabling its use e.g. in other templates.
   *
   * Use this class as in <code>fixed_int_power@<5,2@>::%value</code> to
   * compute 5<sup>2</sup>.
   *
   * @deprecated This template has been deprecated in favor of C++11's support
   * for <code>constexpr</code> calculations, e.g., use
   *
   * @code
   * constexpr int value = Utilities::pow(2, dim);
   * @endcode
   *
   * instead of
   *
   * @code
   * const int value = Utilities::fixed_int_power<2, dim>::value;
   * @endcode
   *
   * to obtain a constant expression for <code>value</code>.
   */
  template <int a, int N>
  struct DEAL_II_DEPRECATED fixed_int_power
  {
    static const int value = a *fixed_int_power<a,N-1>::value;
  };

  /**
   * Base case for the power operation with <code>N=0</code>, which gives the
   * result 1.
   *
   * @deprecated This template is deprecated: see the note in the general
   * version of this template for more information.
   */
  template <int a>
  struct DEAL_II_DEPRECATED fixed_int_power<a,0>
  {
    static const int value = 1;
  };

  /**
   * A replacement for <code>std::pow</code> that allows compile-time
   * calculations for constant expression arguments.
   */
  constexpr
  unsigned int
  pow(const unsigned int base, const unsigned int iexp)
  {
    return iexp == 0 ? 1 : base*dealii::Utilities::pow(base, iexp - 1);
  }

  /**
   * Optimized replacement for <tt>std::lower_bound</tt> for searching within
   * the range of column indices. Slashes execution time by approximately one
   * half for the present application, partly because the binary search is
   * replaced by a linear search for small loop lengths.
   *
   * Another reason for this function is rather obscure: when using the GCC
   * libstdc++ function std::lower_bound, complexity is O(log(N)) as required.
   * However, when using the debug version of the GCC libstdc++ as we do when
   * running the testsuite, then std::lower_bound tests whether the sequence
   * is in fact partitioned with respect to the pivot 'value' (i.e. in essence
   * that the sequence is sorted as required for binary search to work).
   * However, verifying this means that the complexity of std::lower_bound
   * jumps to O(N); we call this function O(N) times below, making the overall
   * complexity O(N**2). The consequence is that a few tests with big meshes
   * completely run off the wall time limit for tests and fail with the
   * libstdc++ debug mode
   *
   * This function simply makes the assumption that the sequence is sorted,
   * and we simply don't do the additional check.
   */
  template <typename Iterator, typename T>
  Iterator
  lower_bound (Iterator  first,
               Iterator  last,
               const T  &val);


  /**
   * The same function as above, but taking an argument that is used to
   * compare individual elements of the sequence of objects pointed to by the
   * iterators.
   */
  template <typename Iterator, typename T, typename Comp>
  Iterator
  lower_bound (Iterator   first,
               Iterator   last,
               const T   &val,
               const Comp comp);

  /**
   * Given a permutation vector (i.e. a vector $p_0\ldots p_{N-1}$ where each
   * $p_i\in [0,N)$ and $p_i\neq p_j$ for $i\neq j$), produce the reverse
   * permutation $q_i=N-1-p_i$.
   */
  std::vector<unsigned int>
  reverse_permutation (const std::vector<unsigned int> &permutation);

  /**
   * Given a permutation vector (i.e. a vector $p_0\ldots p_{N-1}$ where each
   * $p_i\in [0,N)$ and $p_i\neq p_j$ for $i\neq j$), produce the inverse
   * permutation $q_0\ldots q_{N-1}$ so that $q_{p_i}=p_{q_i}=i$.
   */
  std::vector<unsigned int>
  invert_permutation (const std::vector<unsigned int> &permutation);

  /**
   * Given a permutation vector (i.e. a vector $p_0\ldots p_{N-1}$ where each
   * $p_i\in [0,N)$ and $p_i\neq p_j$ for $i\neq j$), produce the reverse
   * permutation $q_i=N-1-p_i$.
   */
  std::vector<unsigned long long int>
  reverse_permutation (const std::vector<unsigned long long int> &permutation);

  /**
   * Given a permutation vector (i.e. a vector $p_0\ldots p_{N-1}$ where each
   * $p_i\in [0,N)$ and $p_i\neq p_j$ for $i\neq j$), produce the inverse
   * permutation $q_0\ldots q_{N-1}$ so that $q_{p_i}=p_{q_i}=i$.
   */
  std::vector<unsigned long long int>
  invert_permutation (const std::vector<unsigned long long int> &permutation);

  /**
   * Given an arbitrary object of type T, use boost::serialization utilities
   * to pack the object into a vector of characters. The object can be unpacked
   * using the Utilities::unpack function below.
   *
   * If the library has been compiled with ZLIB enabled, then the output buffer
   * is compressed.
   *
   * @author Timo Heister, Wolfgang Bangerth, 2017.
   */
  template <typename T>
  std::vector<char> pack (const T &object);

  /**
   * Given a vector of characters, obtained through a call to the function
   * Utilities::pack, restore its content in an object of type T.
   *
   * This function uses boost::serialization utilities to unpack the object
   * from a vector of characters, and it is the inverse of the function
   * Utilities::pack().
   *
   * @note Since no arguments to this function depend on the template type
   *  @p T, you must manually specify the template argument when calling
   *  this function.
   *
   * @note If you want to pack() or unpack() arrays of objects, then the
   *  following works:
   *  @code
   *    double array[3] = {1,2,3};
   *    std::vector<char> buffer = Utilities::pack(array);
   *  @endcode
   *  However, the converse does not:
   *  @code
   *    array = Utilities::unpack<double[3]>(buffer);
   *  @endcode
   *  This is because C++ does not allow functions to return arrays.
   *  Consequently, there is a separate unpack() function for arrays, see
   *  below.
   *
   * @author Timo Heister, Wolfgang Bangerth, 2017.
   */
  template <typename T>
  T unpack (const std::vector<char> &buffer);

  /**
   * Given a vector of characters, obtained through a call to the function
   * Utilities::pack, restore its content in an array of type T.
   *
   * This function uses boost::serialization utilities to unpack the object
   * from a vector of characters, and it is the inverse of the function
   * Utilities::pack().
   *
   * @note This function exists due to a quirk of C++. Specifically,
   *  if you want to pack() or unpack() arrays of objects, then the
   *  following works:
   *  @code
   *    double array[3] = {1,2,3};
   *    std::vector<char> buffer = Utilities::pack(array);
   *  @endcode
   *  However, the converse does not:
   *  @code
   *    array = Utilities::unpack<double[3]>(buffer);
   *  @endcode
   *  This is because C++ does not allow functions to return arrays.
   *  The current function therefore allows to write
   *  @code
   *    Utilities::unpack(buffer, array);
   *  @endcode
   *  Note that unlike the other unpack() function, it is not necessary
   *  to explicitly specify the template arguments since they can be
   *  deduced from the second argument.
   *
   * @author Timo Heister, Wolfgang Bangerth, 2017.
   */
  template <typename T, int N>
  void unpack (const std::vector<char> &buffer,
               T (&unpacked_object)[N]);

  /**
   * A namespace for utility functions that probe system properties.
   *
   * @ingroup utilities
   */
  namespace System
  {

    /**
     * Return the CPU load as returned by "uptime". Note that the
     * interpretation of this number depends on the actual number of
     * processors in the machine. This is presently only implemented on Linux,
     * using the /proc/loadavg pseudo-file, on other systems we simply return
     * zero.
     */
    double get_cpu_load ();

    /**
     * Structure that hold information about memory usage in kB. Used by
     * get_memory_stats(). See man 5 proc entry /status for details.
     */
    struct MemoryStats
    {
      /**
       * Peak virtual memory size in kB.
       */
      unsigned long int VmPeak;

      /**
       * Current virtual memory size in kB.
       */
      unsigned long int VmSize;

      /**
       * Peak resident memory size in kB. Also known as "high water mark" (HWM).
       */
      unsigned long int VmHWM;

      /**
       * Current resident memory size in kB. Also known as "resident set size" (RSS).
       */
      unsigned long int VmRSS;
    };


    /**
     * Fill the @p stats structure with information about the memory
     * consumption of this process. This is only implemented on Linux.
     */
    void get_memory_stats (MemoryStats &stats);


    /**
     * Return the name of the host this process runs on.
     */
    std::string get_hostname ();


    /**
     * Return the present time as HH:MM:SS.
     */
    std::string get_time ();

    /**
     * Return the present date as YYYY/MM/DD. MM and DD may be either one or
     * two digits.
     */
    std::string get_date ();

    /**
     * Call the system function posix_memalign, or a replacement function if
     * not available, to allocate memory with a certain minimal alignment. The
     * first argument will then return a pointer to this memory block that can
     * be released later on through a standard <code>free</code> call.
     *
     * @param memptr The address of a pointer variable that will after this
     * call point to the allocated memory.
     * @param alignment The minimal alignment of the memory block, in bytes.
     * @param size The size of the memory block to be allocated, in bytes.
     *
     * @note This function checks internally for error codes, rather than
     * leaving this task to the calling site.
     */
    void posix_memalign (void **memptr, size_t alignment, size_t size);
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
     * Return a Trilinos Epetra_Comm object needed for creation of
     * Epetra_Maps.
     *
     * If deal.II has been configured to use a compiler that does not support
     * MPI then the resulting communicator will be a serial one. Otherwise,
     * the communicator will correspond to MPI_COMM_WORLD, i.e. a communicator
     * that encompasses all processes within this MPI universe.
     */
    const Epetra_Comm &comm_world();

    /**
     * Return a Trilinos Epetra_Comm object needed for creation of
     * Epetra_Maps.
     *
     * If deal.II has been configured to use a compiler that does not support
     * MPI then the resulting communicator will be a serial one. Otherwise,
     * the communicator will correspond to MPI_COMM_SELF, i.e. a communicator
     * that comprises only this one processor.
     */
    const Epetra_Comm &comm_self();

    /**
     * Given a communicator, duplicate it. If the given communicator is
     * serial, that means to just return a copy of itself. On the other hand,
     * if it is %parallel, we duplicate the underlying MPI_Comm object: we
     * create a separate MPI communicator that contains the same processors
     * and in the same order but has a separate identifier distinct from the
     * given communicator. The function returns a pointer to a new object of a
     * class derived from Epetra_Comm. The caller of this function needs to
     * assume ownership of this function. The returned object should be
     * destroyed using the destroy_communicator() function.
     *
     * This facility is used to separate streams of communication. For
     * example, a program could simply use MPI_Comm_World for everything. But
     * it is easy to come up with scenarios where sometimes not all processors
     * participate in a communication that is intended to be global -- for
     * example if we assemble a matrix on a coarse mesh with fewer cells than
     * there are processors, some processors may not sync their matrices with
     * the rest because they haven't written into it because they own no
     * cells. That's clearly a bug. However, if these processors just continue
     * their work, and the next %parallel operation happens to be a sync on a
     * different matrix, then the sync could succeed -- by accident, since
     * different processors are talking about different matrices.
     *
     * This kind of situation can be avoided if we use different communicators
     * for different matrices which reduces the likelihood that communications
     * meant to be separate aren't recognized as such just because they happen
     * on the same communicator. In addition, it is conceivable that some MPI
     * operations can be parallelized using multiple threads because their
     * communicators identifies the communication in question, not their
     * relative timing as is the case in a sequential program that just uses a
     * single communicator.
     */
    Epetra_Comm *
    duplicate_communicator (const Epetra_Comm &communicator);

    /**
     * Given an Epetra communicator that was created by the
     * duplicate_communicator() function, destroy the underlying MPI
     * communicator object and reset the Epetra_Comm object to a the result of
     * comm_self().
     *
     * It is necessary to call this function at the time when the result of
     * duplicate_communicator() is no longer needed. The reason is that in
     * that function, we first create a new MPI_Comm object and then create an
     * Epetra_Comm around it. While we can take care of destroying the latter,
     * it doesn't destroy the communicator since it can only assume that it
     * may also be still used by other objects in the program. Consequently,
     * we have to take care of destroying it ourselves, explicitly.
     *
     * This function does exactly that. Because this has to happen while the
     * Epetra_Comm object is still around, it first resets the latter and then
     * destroys the communicator object.
     *
     * @note If you call this function on an Epetra_Comm object that is not
     * created by duplicate_communicator(), you are likely doing something
     * quite wrong. Don't do this.
     */
    void
    destroy_communicator (Epetra_Comm &communicator);

    /**
     * Return the number of MPI processes there exist in the given
     * @ref GlossMPICommunicator "communicator"
     * object. If this is a sequential job (i.e., the program
     * is not using MPI at all, or is using MPI but has been started with
     * only one MPI process), then the communicator necessarily involves
     * only one process and the function returns 1.
     */
    unsigned int get_n_mpi_processes (const Epetra_Comm &mpi_communicator);

    /**
     * Return the number of the present MPI process in the space of processes
     * described by the given communicator. This will be a unique value for
     * each process between zero and (less than) the number of all processes
     * (given by get_n_mpi_processes()).
     */
    unsigned int get_this_mpi_process (const Epetra_Comm &mpi_communicator);

    /**
     * Given a Trilinos Epetra map, create a new map that has the same
     * subdivision of elements to processors but uses the given communicator
     * object instead of the one stored in the first argument. In essence,
     * this means that we create a map that communicates among the same
     * processors in the same way, but using a separate channel.
     *
     * This function is typically used with a communicator that has been
     * obtained by the duplicate_communicator() function.
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
    Assert (N>=0, ExcNotImplemented());
    switch (N)
      {
      case 0:
        return dealii::internal::NumberType<T>::value(1);
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



  template <typename Iterator, typename T>
  inline
  Iterator
  lower_bound (Iterator  first,
               Iterator  last,
               const T  &val)
  {
    return Utilities::lower_bound (first, last, val,
                                   std::less<T>());
  }



  template <typename Iterator, typename T, typename Comp>
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
                DEAL_II_FALLTHROUGH;
              case 6:
                if (!comp(*first, val))
                  return first;
                ++first;
                DEAL_II_FALLTHROUGH;
              case 5:
                if (!comp(*first, val))
                  return first;
                ++first;
                DEAL_II_FALLTHROUGH;
              case 4:
                if (!comp(*first, val))
                  return first;
                ++first;
                DEAL_II_FALLTHROUGH;
              case 3:
                if (!comp(*first, val))
                  return first;
                ++first;
                DEAL_II_FALLTHROUGH;
              case 2:
                if (!comp(*first, val))
                  return first;
                ++first;
                DEAL_II_FALLTHROUGH;
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
                // shouldn't have
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



  template <typename T>
  std::vector<char> pack(const T &object)
  {
    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
    //
    // we have to work around the fact that GCC 4.8.x claims to be C++
    // conforming, but is not actually as it does not implement
    // std::is_trivially_copyable.
    if (
#if __GNUG__ && __GNUC__ < 5
      __has_trivial_copy(T)
#else
      std::is_trivially_copyable<T>()
#endif
      &&
      sizeof(T)<256)
      {
        std::vector<char> buffer (sizeof(T));
        std::memcpy (buffer.data(), &object, sizeof(T));
        return buffer;
      }
    else
      {
        // set up a buffer and then use it as the target of a compressing
        // stream into which we serialize the current object
        std::vector<char> buffer;
        {
#ifdef DEAL_II_WITH_ZLIB
          boost::iostreams::filtering_ostream out;
          out.push(boost::iostreams::gzip_compressor
                   (boost::iostreams::gzip_params
                    (boost::iostreams::gzip::best_compression)));
          out.push(boost::iostreams::back_inserter(buffer));

          boost::archive::binary_oarchive archive(out);
          archive << object;
          out.flush();
#else
          std::ostringstream out;
          boost::archive::binary_oarchive archive(out);
          archive << object;
          const std::string &s = out.str();
          buffer.reserve(s.size());
          buffer.assign(s.begin(), s.end());
#endif
        }
        return buffer;
      }
  }


  template <typename T>
  T unpack(const std::vector<char> &buffer)
  {
    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
    //
    // we have to work around the fact that GCC 4.8.x claims to be C++
    // conforming, but is not actually as it does not implement
    // std::is_trivially_copyable.
    if (
#if __GNUG__ && __GNUC__ < 5
      __has_trivial_copy(T)
#else
      std::is_trivially_copyable<T>()
#endif
      &&
      sizeof(T)<256)
      {
        Assert (buffer.size() == sizeof(T), ExcInternalError());
        T object;
        std::memcpy (&object, buffer.data(), sizeof(T));
        return object;
      }
    else
      {
        std::string decompressed_buffer;
        T object;

        // first decompress the buffer
        {
#ifdef DEAL_II_WITH_ZLIB
          boost::iostreams::filtering_ostream decompressing_stream;
          decompressing_stream.push(boost::iostreams::gzip_decompressor());
          decompressing_stream.push(boost::iostreams::back_inserter(decompressed_buffer));
          decompressing_stream.write (buffer.data(), buffer.size());
#else
          decompressed_buffer.assign (buffer.begin(), buffer.end());
#endif
        }

        // then restore the object from the buffer
        std::istringstream in(decompressed_buffer);
        boost::archive::binary_iarchive archive(in);

        archive >> object;
        return object;
      }
  }



  template <typename T, int N>
  void unpack (const std::vector<char> &buffer,
               T (&unpacked_object)[N])
  {
    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
    //
    // we have to work around the fact that GCC 4.8.x claims to be C++
    // conforming, but is not actually as it does not implement
    // std::is_trivially_copyable.
    if (
#if __GNUG__ && __GNUC__ < 5
      __has_trivial_copy(T)
#else
      std::is_trivially_copyable<T>()
#endif
      &&
      sizeof(T)*N<256)
      {
        Assert (buffer.size() == sizeof(T)*N, ExcInternalError());
        std::memcpy (unpacked_object, buffer.data(), sizeof(T)*N);
      }
    else
      {
        std::string decompressed_buffer;

        // first decompress the buffer
        {
#ifdef DEAL_II_WITH_ZLIB
          boost::iostreams::filtering_ostream decompressing_stream;
          decompressing_stream.push(boost::iostreams::gzip_decompressor());
          decompressing_stream.push(boost::iostreams::back_inserter(decompressed_buffer));
          decompressing_stream.write (buffer.data(), buffer.size());
#else
          decompressed_buffer.assign (buffer.begin(), buffer.end());
#endif
        }

        // then restore the object from the buffer
        std::istringstream in(decompressed_buffer);
        boost::archive::binary_iarchive archive(in);

        archive >> unpacked_object;
      }
  }

}


DEAL_II_NAMESPACE_CLOSE

#ifndef DOXYGEN
namespace boost
{
  namespace serialization
  {

    // Provides boost and c++11 with a way to serialize tuples and pairs automatically
    template<int N>
    struct Serialize
    {
      template<class Archive, typename... Args>
      static void serialize(Archive &ar, std::tuple<Args...> &t, const unsigned int version)
      {
        ar &std::get<N-1>(t);
        Serialize<N-1>::serialize(ar, t, version);
      }
    };

    template<>
    struct Serialize<0>
    {
      template<class Archive, typename... Args>
      static void serialize(Archive &ar, std::tuple<Args...> &t, const unsigned int version)
      {
        (void) ar;
        (void) t;
        (void) version;
      }
    };

    template<class Archive, typename... Args>
    void serialize(Archive &ar, std::tuple<Args...> &t, const unsigned int version)
    {
      Serialize<sizeof...(Args)>::serialize(ar, t, version);
    }
  }
}
#endif

#endif
