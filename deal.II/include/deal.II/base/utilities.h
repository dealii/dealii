//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__utilities_h
#define __deal2__utilities_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <vector>
#include <utility>
#include <functional>
#include <string>

#if defined(DEAL_II_COMPILER_SUPPORTS_MPI) || defined(DEAL_II_USE_PETSC)
#include <mpi.h>
#else

extern boost::mutex m;typedef int MPI_Comm;
#endif

#ifdef DEAL_II_USE_TRILINOS
#  include <Epetra_Comm.h>
#  include <Epetra_Map.h>
#  ifdef DEAL_II_COMPILER_SUPPORTS_MPI
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
                                    * positions marked by @p delimiter in the text.
                                    * If this is not possible, return the shortest
                                    * lines that are longer than @p width.
                                    * The default value of the delimiter is a
                                    * space character.
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
                                    * A namespace for utility functions that
                                    * abstract certain operations using the
				    * Message Passing Interface (MPI) or
				    * provide fallback operations in
				    * case deal.II is configured not to use
				    * MPI at all.
                                    *
                                    * @ingroup utilities
                                    */
  namespace MPI
  {
				     /**
				      * Return the number of MPI processes
				      * there exist in the given communicator
				      * object. If this is a sequential job,
				      * it returns 1.
				      */
    unsigned int n_mpi_processes (const MPI_Comm &mpi_communicator);

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
    unsigned int this_mpi_process (const MPI_Comm &mpi_communicator);

				     /**
				      * Consider an unstructured
				      * communication pattern where
				      * every process in an MPI
				      * universe wants to send some
				      * data to a subset of the other
				      * processors. To do that, the
				      * other processors need to know
				      * who to expect messages
				      * from. This function computes
				      * this information.
				      *
				      * @param mpi_comm A communicator
				      * that describes the processors
				      * that are going to communicate
				      * with each other.
				      *
				      * @param destinations The list
				      * of processors the current
				      * process wants to send
				      * information to. This list need
				      * not be sorted in any way. If
				      * it contains duplicate entries
				      * that means that multiple
				      * messages are intended for a
				      * given destination.
				      *
				      * @return A list of processors
				      * that have indicated that they
				      * want to send something to the
				      * current processor. The
				      * resulting list is not
				      * sorted. It may contain
				      * duplicate entries if
				      * processors enter the same
				      * destination more than once in
				      * their destinations list.
				      */
    std::vector<unsigned int>
    compute_point_to_point_communication_pattern (const MPI_Comm & mpi_comm,
						  const std::vector<unsigned int> & destinations);

				     /**
				      * Given a communicator, generate a new
				      * communicator that contains the same
				      * set of processors but that has a
				      * different, unique identifier.
				      *
				      * This functionality can be used to
				      * ensure that different objects, such as
				      * distributed matrices, each have unique
				      * communicators over which they can
				      * interact without interfering with each
				      * other.
				      *
				      * When no longer needed, the
				      * communicator created here needs to
				      * be destroyed using
				      * <code>MPI_Comm_free</code>.
				      */
    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator);

    /**
     * Return the sum over all processors of the value @p t. This function
     * is collective over all processors given in the communicator. If
     * deal.II is not configured for use of MPI, this function simply
     * returns the value of @p t. This function corresponds to the
     * <code>MPI_Allreduce</code> function, i.e. all processors receive
     * the result of this operation.
     *
     * @note This function is only implemented for certain template
     * arguments <code>T</code>, namely <code>float, double, int,
     * unsigned int</code>.
     */
    template <typename T>
    T sum (const T &t,
	   const MPI_Comm &mpi_communicator);

				     /**
				      * Data structure to store the result of
				      * calculate_collective_mpi_min_max_avg.
				      */
    struct MinMaxAvg
    {
	double sum;
	double min;
	double max;
	unsigned int min_index;
	unsigned int max_index;
	double avg;
    };

				     /**
				      * Returns sum, average, minimum,
				      * maximum, processor id of minimum and
				      * maximum as a collective operation of
				      * on the given MPI communicator @param
				      * mpi_communicator . Each processor's
				      * value is given in @param my_value and
				      * the result will be returned in @param
				      * result . The result is available on all
				      * machines.
				      */
    MinMaxAvg
    min_max_avg (const double my_value,
		 const MPI_Comm &mpi_communicator);
  }
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
				      * Fills the @param stats structure with
				      * information about the memory
				      * consumption of this process. This is
				      * only implemented on Linux.
				      */
    void get_memory_stats (MemoryStats & stats);


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
    bool program_uses_mpi();

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
    unsigned int get_n_mpi_processes (const MPI_Comm &mpi_communicator);

				     /**
				      * This function is an alias for
				      * Utilities::MPI::this_mpi_process.
				      *
				      * @deprecated
				      */
    unsigned int get_this_mpi_process (const MPI_Comm &mpi_communicator);


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
					  MinMaxAvg & result);

				     /** @} */


				     /**
				      * A class that is used to initialize the
				      * MPI system at the beginning of a
				      * program and to shut it down again at
				      * the end.
				      *
				      * If a program uses MPI one would
				      * typically just create an object of
				      * this type at the beginning of
				      * <code>main()</code>. The constructor
				      * of this class then runs
				      * <code>MPI_Init()</code> with the given
				      * arguments. At the end of the program,
				      * the compiler will invoke the
				      * destructor of this object which in
				      * turns calls <code>MPI_Finalize</code>
				      * to shut down the MPI system.
				      *
				      * This class is used in @ref step_32
				      * "step-32", for example.
				      */
    class MPI_InitFinalize
    {
      public:
					 /**
					  * Constructor. Takes the arguments
					  * from the command line (in case of
					  * MPI, the number of processes is
					  * specified there), and sets up a
					  * respective communicator by calling
					  * <tt>MPI_Init()</tt>. This
					  * constructor can only be called once
					  * in a program, since MPI cannot be
					  * initialized twice.
					  */
	MPI_InitFinalize (int    &argc,
			  char** &argv);

					 /**
					  * Destructor. Calls
					  * <tt>MPI_Finalize()</tt> in
					  * case this class owns the MPI
					  * process.
					  */
	~MPI_InitFinalize();

      private:
					 /**
					  * This flag tells the class
					  * whether it owns the MPI
					  * process (i.e., it has been
					  * constructed using the
					  * argc/argv input, or it has
					  * been copied). In the former
					  * case, the command
					  * <tt>MPI_Finalize()</tt> will
					  * be called at destruction.
					  */
	const bool owns_mpi;
    };
  }


  /**
   * This class provides the basic structures for the use of the Trilinos
   * classes such as matrices, vectors, and preconditioners. The most
   * important function in this class is <tt>comm()</tt>, which is needed
   * for the initialization of Trilinos Epetra_Maps, which design the
   * %parallel distribution of vectors and matrices. Moreover, this class
   * provides a unified interface to both serial and %parallel
   * implementations of Trilinos, sets up the MPI communicator in case the
   * programs are run in %parallel, and correctly terminates all processes
   * when the destructor is called. An example usage of this class is shown
   * in the tutorial program step-32.
   */
#ifdef DEAL_II_USE_TRILINOS
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
    const Epetra_Comm& comm_world();

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
    const Epetra_Comm& comm_self();

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
        case 0:
	      return 1;
	case 1:
	      return n;
	case 2:
	      return n*n;
	case 3:
	      return n*n*n;
	case 4:
	      return n*n*n*n;
	default:
	      T result = 1;
	      for (int d=0;d<N;++d)
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
    unsigned int len = last-first;

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


  namespace MPI
  {
    namespace internal
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      /**
       * Return the corresponding MPI data type id for the argument given.
       */
      inline MPI_Datatype mpi_type_id (const int *)
      {
	return MPI_INT;
      }


      inline MPI_Datatype mpi_type_id (const unsigned int *)
      {
	return MPI_UNSIGNED;
      }


      inline MPI_Datatype mpi_type_id (const float *)
      {
	return MPI_FLOAT;
      }


      inline MPI_Datatype mpi_type_id (const double *)
      {
	return MPI_DOUBLE;
      }
#endif
    }

    template <typename T>
    inline
    T sum (const T &t,
	   const MPI_Comm &mpi_communicator)
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      T sum;
      MPI_Allreduce (const_cast<void*>(static_cast<const void*>(&t)),
		     &sum, 1, internal::mpi_type_id(&t), MPI_SUM,
		     mpi_communicator);
      return sum;
#else
      (void)mpi_communicator;
      return t;
#endif
    }
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif

