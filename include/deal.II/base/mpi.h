// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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

#ifndef __deal2__mpi_h
#define __deal2__mpi_h

#include <deal.II/base/config.h>
#include <vector>

#if defined(DEAL_II_WITH_MPI) || defined(DEAL_II_WITH_PETSC)
#  include <mpi.h>
// Check whether <mpi.h> is a suitable
// include for us (if MPI_SEEK_SET is not
// defined, we'll die anyway):
#  ifndef MPI_SEEK_SET
#    error "The buildsystem included an insufficient mpi.h header that does not export MPI_SEEK_SET"
#  endif

#else
// without MPI, we would still like to use
// some constructs with MPI data
// types. Therefore, create some dummies
typedef int MPI_Comm;
const int MPI_COMM_SELF = 0;
typedef int MPI_Datatype;
typedef int MPI_Op;
namespace MPI
{
  static const unsigned int UNSIGNED = 0;
  static const unsigned int LONG_DOUBLE = 0;
  static const unsigned int LONG_DOUBLE_COMPLEX = 0;
  static const unsigned int MAX = 0;
  static const unsigned int MIN = 0;
  static const unsigned int SUM = 0;
}
#endif

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  /**
   * A namespace for utility functions that abstract certain operations using
   * the Message Passing Interface (MPI) or provide fallback operations in
   * case deal.II is configured not to use MPI at all.
   *
   * @ingroup utilities
   */
  namespace MPI
  {
    /**
     * Return the number of MPI processes there exist in the given
     * communicator object. If this is a sequential job, it returns 1.
     */
    unsigned int n_mpi_processes (const MPI_Comm &mpi_communicator);

    /**
     * Return the number of the present MPI process in the space of processes
     * described by the given communicator. This will be a unique value for
     * each process between zero and (less than) the number of all processes
     * (given by get_n_mpi_processes()).
     */
    unsigned int this_mpi_process (const MPI_Comm &mpi_communicator);

    /**
     * Consider an unstructured communication pattern where every process in
     * an MPI universe wants to send some data to a subset of the other
     * processors. To do that, the other processors need to know who to expect
     * messages from. This function computes this information.
     *
     * @param mpi_comm A communicator that describes the processors that are
     * going to communicate with each other.
     *
     * @param destinations The list of processors the current process wants to
     * send information to. This list need not be sorted in any way. If it
     * contains duplicate entries that means that multiple messages are
     * intended for a given destination.
     *
     * @return A list of processors that have indicated that they want to send
     * something to the current processor. The resulting list is not sorted.
     * It may contain duplicate entries if processors enter the same
     * destination more than once in their destinations list.
     */
    std::vector<unsigned int>
    compute_point_to_point_communication_pattern (const MPI_Comm &mpi_comm,
                                                  const std::vector<unsigned int> &destinations);

    /**
     * Given a communicator, generate a new communicator that contains the
     * same set of processors but that has a different, unique identifier.
     *
     * This functionality can be used to ensure that different objects, such
     * as distributed matrices, each have unique communicators over which they
     * can interact without interfering with each other.
     *
     * When no longer needed, the communicator created here needs to be
     * destroyed using <code>MPI_Comm_free</code>.
     */
    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator);

    /**
     * Return the sum over all processors of the value @p t. This function is
     * collective over all processors given in the communicator. If deal.II is
     * not configured for use of MPI, this function simply returns the value
     * of @p t. This function corresponds to the <code>MPI_Allreduce</code>
     * function, i.e. all processors receive the result of this operation.
     *
     * @note Sometimes, not all processors need a results and in that case one
     * would call the <code>MPI_Reduce</code> function instead of the
     * <code>MPI_Allreduce</code> function. The latter is at most twice as
     * expensive, so if you are concerned about performance, it may be
     * worthwhile investigating whether your algorithm indeed needs the result
     * everywhere or whether you could get away with calling the current
     * function and getting the result everywhere.
     *
     * @note This function is only implemented for certain template arguments
     * <code>T</code>, namely <code>float, double, int, unsigned int</code>.
     */
    template <typename T>
    T sum (const T &t,
           const MPI_Comm &mpi_communicator);

    /**
     * Like the previous function, but take the sums over the elements of an
     * array of length N. In other words, the i-th element of the results
     * array is the sum over the i-th entries of the input arrays from each
     * processor.
     */
    template <typename T, unsigned int N>
    inline
    void sum (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&sums)[N]);

    /**
     * Like the previous function, but take the sums over the elements of a
     * std::vector. In other words, the i-th element of the results array is
     * the sum over the i-th entries of the input arrays from each processor.
     */
    template <typename T>
    inline
    void sum (const std::vector<T> &values,
              const MPI_Comm &mpi_communicator,
              std::vector<T> &sums);

    /**
     * Return the maximum over all processors of the value @p t. This function
     * is collective over all processors given in the communicator. If deal.II
     * is not configured for use of MPI, this function simply returns the
     * value of @p t. This function corresponds to the
     * <code>MPI_Allreduce</code> function, i.e. all processors receive the
     * result of this operation.
     *
     * @note Sometimes, not all processors need a results and in that case one
     * would call the <code>MPI_Reduce</code> function instead of the
     * <code>MPI_Allreduce</code> function. The latter is at most twice as
     * expensive, so if you are concerned about performance, it may be
     * worthwhile investigating whether your algorithm indeed needs the result
     * everywhere or whether you could get away with calling the current
     * function and getting the result everywhere.
     *
     * @note This function is only implemented for certain template arguments
     * <code>T</code>, namely <code>float, double, int, unsigned int</code>.
     */
    template <typename T>
    T max (const T &t,
           const MPI_Comm &mpi_communicator);

    /**
     * Like the previous function, but take the maxima over the elements of an
     * array of length N. In other words, the i-th element of the results
     * array is the maximum of the i-th entries of the input arrays from each
     * processor.
     */
    template <typename T, unsigned int N>
    inline
    void max (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&maxima)[N]);

    /**
     * Like the previous function, but take the maximum over the elements of a
     * std::vector. In other words, the i-th element of the results array is
     * the maximum over the i-th entries of the input arrays from each
     * processor.
     */
    template <typename T>
    inline
    void max (const std::vector<T> &values,
              const MPI_Comm &mpi_communicator,
              std::vector<T> &maxima);

    /**
     * Data structure to store the result of min_max_avg().
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
     * Returns sum, average, minimum, maximum, processor id of minimum and
     * maximum as a collective operation of on the given MPI communicator @p
     * mpi_communicator . Each processor's value is given in @p my_value and
     * the result will be returned. The result is available on all machines.
     *
     * @note Sometimes, not all processors need a results and in that case one
     * would call the <code>MPI_Reduce</code> function instead of the
     * <code>MPI_Allreduce</code> function. The latter is at most twice as
     * expensive, so if you are concerned about performance, it may be
     * worthwhile investigating whether your algorithm indeed needs the result
     * everywhere or whether you could get away with calling the current
     * function and getting the result everywhere.
     */
    MinMaxAvg
    min_max_avg (const double my_value,
                 const MPI_Comm &mpi_communicator);



    /**
     * A class that is used to initialize the MPI system at the beginning of a
     * program and to shut it down again at the end. It also allows you to
     * control the number threads used in each MPI task.
     *
     * If deal.II is configured with PETSc, the library will also be
     * initialized in the beginning and destructed at the end automatically
     * (internally by calling PetscInitialize() and PetscFinalize()).
     *
     * If a program uses MPI one would typically just create an object of this
     * type at the beginning of <code>main()</code>. The constructor of this
     * class then runs <code>MPI_Init()</code> with the given arguments. At
     * the end of the program, the compiler will invoke the destructor of this
     * object which in turns calls <code>MPI_Finalize</code> to shut down the
     * MPI system.
     *
     * This class is used in step-32, for example.
     */
    class MPI_InitFinalize
    {
    public:

      /**
       * Constructor. Takes the arguments from the command line (in case of
       * MPI, the number of processes is specified there), and sets up a
       * respective communicator by calling <tt>MPI_Init()</tt>. This
       * constructor can only be called once in a program, since MPI cannot be
       * initialized twice.
       *
       * This constructor sets max_num_threads to 1 (see other constructor).
       */
      MPI_InitFinalize (int    &argc,
                        char ** &argv) /*DEAL_II_DEPRECATED*/;

      /**
       * Initialize MPI (and, if deal.II was configured to use it, PETSc)
       * and set the number of threads used by deal.II (via the underlying
       * Threading Building Blocks library) to the given parameter.
       *
       * @param[in,out] argc A reference to the 'argc' argument passed to main. This
       *   argument is used to initialize MPI (and, possibly, PETSc) as they
       *   read arguments from the command line.
       * @param[in,out] argv A reference to the 'argv' argument passed to main.
       * @param[in] max_num_threads The maximal number of threads this MPI process
       *   should utilize. If this argument is set to
       *   numbers::invalid_unsigned_int, the number of threads is determined by
       *   automatically in the following way: the number of
       *   threads to run on this MPI process is set in such a way that all of
       *   the cores in your node are spoken for. In other words, if you
       *   have started one MPI process per node, setting this argument is
       *   equivalent to setting it to the number of cores present in the node
       *   this MPI process runs on. If you have started as many MPI
       *   processes per node as there are cores on each node, then
       *   this is equivalent to passing 1 as the argument. On the
       *   other hand, if, for example, you start 4 MPI processes
       *   on each 16-core node, then this option will start 4 worker
       *   threads for each node. If you start 3 processes on an 8 core
       *   node, then they will start 3, 3 and 2 threads, respectively.
       *
       * @note This function calls MultithreadInfo::set_thread_limit()
       * with either @p max_num_threads or, following the discussion above, a
       * number of threads equal to the number of cores allocated to this
       * MPI process. However, MultithreadInfo::set_thread_limit() in turn also
       * evaluates the environment variable DEAL_II_NUM_THREADS. Finally, the worker
       * threads can only be created on cores to which the current MPI process has
       * access to; some MPI implementations limit the number of cores each process
       * has access to to one or a subset of cores in order to ensure better cache
       * behavior. Consequently, the number of threads that will really be created
       * will be the minimum of the argument passed here, the environment variable
       * (if set), and the number of cores accessible to the thread.
       */
      MPI_InitFinalize (int    &argc,
                        char ** &argv,
                        const unsigned int max_num_threads);

      /**
       * Destructor. Calls <tt>MPI_Finalize()</tt> in case this class owns the
       * MPI process.
       */
      ~MPI_InitFinalize();

    private:
      /**
       * This flag tells the class whether it owns the MPI process (i.e., it
       * has been constructed using the argc/argv input, or it has been
       * copied). In the former case, the command <tt>MPI_Finalize()</tt> will
       * be called at destruction.
       */
      const bool owns_mpi;


      /**
       * A common function called by all of the constructors.
       */
      void do_init(int    &argc,
                   char ** &argv);
    };

    namespace internal
    {
#ifdef DEAL_II_WITH_MPI
      /**
       * Return the corresponding MPI data type id for the argument given.
       */
      inline MPI_Datatype mpi_type_id (const int *)
      {
        return MPI_INT;
      }


      inline MPI_Datatype mpi_type_id (const long int *)
      {
        return MPI_LONG;
      }


      inline MPI_Datatype mpi_type_id (const unsigned int *)
      {
        return MPI_UNSIGNED;
      }


      inline MPI_Datatype mpi_type_id (const unsigned long int *)
      {
        return MPI_UNSIGNED_LONG;
      }


      inline MPI_Datatype mpi_type_id (const unsigned long long int *)
      {
        return MPI_UNSIGNED_LONG_LONG;
      }


      inline MPI_Datatype mpi_type_id (const float *)
      {
        return MPI_FLOAT;
      }


      inline MPI_Datatype mpi_type_id (const double *)
      {
        return MPI_DOUBLE;
      }


      inline MPI_Datatype mpi_type_id (const long double *)
      {
        return MPI_LONG_DOUBLE;
      }
#endif
    }


    template <typename T>
    inline
    T sum (const T &t,
           const MPI_Comm &mpi_communicator)
    {
#ifdef DEAL_II_WITH_MPI
      T sum;
      MPI_Allreduce (const_cast<void *>(static_cast<const void *>(&t)),
                     &sum, 1, internal::mpi_type_id(&t), MPI_SUM,
                     mpi_communicator);
      return sum;
#else
      (void)mpi_communicator;
      return t;
#endif
    }


    template <typename T, unsigned int N>
    inline
    void sum (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&sums)[N])
    {
#ifdef DEAL_II_WITH_MPI
      MPI_Allreduce (const_cast<void *>(static_cast<const void *>(&values[0])),
                     &sums[0], N, internal::mpi_type_id(values), MPI_SUM,
                     mpi_communicator);
#else
      (void)mpi_communicator;
      for (unsigned int i=0; i<N; ++i)
        sums[i] = values[i];
#endif
    }


    template <typename T>
    inline
    void sum (const std::vector<T> &values,
              const MPI_Comm       &mpi_communicator,
              std::vector<T>       &sums)
    {
#ifdef DEAL_II_WITH_MPI
      sums.resize (values.size());
      MPI_Allreduce (const_cast<void *>(static_cast<const void *>(&values[0])),
                     &sums[0], values.size(), internal::mpi_type_id((T *)0), MPI_SUM,
                     mpi_communicator);
#else
      (void)mpi_communicator;
      sums = values;
#endif
    }


    template <typename T>
    inline
    T max (const T &t,
           const MPI_Comm &mpi_communicator)
    {
#ifdef DEAL_II_WITH_MPI
      T sum;
      MPI_Allreduce (const_cast<void *>(static_cast<const void *>(&t)),
                     &sum, 1, internal::mpi_type_id(&t), MPI_MAX,
                     mpi_communicator);
      return sum;
#else
      (void)mpi_communicator;
      return t;
#endif
    }


    template <typename T, unsigned int N>
    inline
    void max (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&maxima)[N])
    {
#ifdef DEAL_II_WITH_MPI
      MPI_Allreduce (const_cast<void *>(static_cast<const void *>(&values[0])),
                     &maxima[0], N, internal::mpi_type_id(values), MPI_MAX,
                     mpi_communicator);
#else
      (void)mpi_communicator;
      for (unsigned int i=0; i<N; ++i)
        maxima[i] = values[i];
#endif
    }


    template <typename T>
    inline
    void max (const std::vector<T> &values,
              const MPI_Comm       &mpi_communicator,
              std::vector<T>       &maxima)
    {
#ifdef DEAL_II_WITH_MPI
      maxima.resize (values.size());
      MPI_Allreduce (const_cast<void *>(static_cast<const void *>(&values[0])),
                     &maxima[0], values.size(), internal::mpi_type_id((T *)0), MPI_MAX,
                     mpi_communicator);
#else
      (void)mpi_communicator;
      maxima = values;
#endif
    }
  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
