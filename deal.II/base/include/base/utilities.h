//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__utilities_h
#define __deal2__utilities_h

#include <base/config.h>
#include <base/exceptions.h>

#include <vector>
#include <utility>
#include <string>

#if defined(DEAL_II_COMPILER_SUPPORTS_MPI) || defined(DEAL_II_USE_PETSC)
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#ifdef DEAL_II_USE_TRILINOS
#  include <Teuchos_RCP.hpp>
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
                                    * lines than are longer than @p width.
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
                                      * Return the name of the host this
                                      * process runs on.
                                      */
    std::string get_hostname ();


                                     /**
                                      * Return the present time as HH:MM:SS.
                                      */
    std::string get_time ();

				     /**
				      * Return the number of MPI processes
				      * there exist in the given communicator
				      * object. If this is a sequential job,
				      * it returns 1.
				      */
    unsigned int get_n_mpi_processes (const MPI_Comm &mpi_communicator);

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
    unsigned int get_this_mpi_process (const MPI_Comm &mpi_communicator);

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
				      */
    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator);

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

				     /**
				      * Returns whether deal.II has been
				      * configured to use MPI and if so
				      * whether MPI has already been
				      * initialized using MPI_Init().
				      */
    bool program_uses_mpi();
  }


  /**
   * This class provides the basic structures for the use of the Trilinos
   * classes such as matrices, vectors, and preconditioners. The most
   * important function in this class is <tt>comm()</tt>, which is needed
   * for the initialization of Trilinos Epetra_Maps, which design the
   * parallel distribution of vectors and matrices. Moreover, this class
   * provides a unified interface to both serial and parallel
   * implementations of Trilinos, sets up the MPI communicator in case the
   * programs are run in parallel, and correctly terminates all processes
   * when the destructor is called. An example usage of this class is shown
   * in the tutorial program @ref step_32 "step-32".
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
				      * parallel, we duplicate the underlying
				      * MPI_Comm object: we create a separate
				      * MPI communicator that contains the
				      * same processors and in the same order
				      * but has a separate identifier distinct
				      * from the given communicator. The
				      * function returns a pointer to a new
				      * object of a class derived from
				      * Epetra_Comm. The caller of this
				      * function needs to assume ownership of
				      * this function.
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
				      * and the next parallel operation
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
    duplicate_map (const Epetra_Map  &map,
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
	      for (unsigned int d=1;d<N;++d)
		result *= n;
	      return result;
      }
  }
  
}


DEAL_II_NAMESPACE_CLOSE

#endif

