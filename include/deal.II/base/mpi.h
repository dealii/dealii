// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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

#ifndef dealii_mpi_h
#define dealii_mpi_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi_tags.h>
#include <deal.II/base/numbers.h>

#include <boost/signals2.hpp>

#include <map>
#include <numeric>
#include <set>
#include <vector>

#if !defined(DEAL_II_WITH_MPI) && !defined(DEAL_II_WITH_PETSC)
// without MPI, we would still like to use
// some constructs with MPI data
// types. Therefore, create some dummies
using MPI_Comm     = int;
using MPI_Request  = int;
using MPI_Datatype = int;
using MPI_Op       = int;
#  ifndef MPI_COMM_WORLD
#    define MPI_COMM_WORLD 0
#  endif
#  ifndef MPI_COMM_SELF
#    define MPI_COMM_SELF 0
#  endif
#  ifndef MPI_REQUEST_NULL
#    define MPI_REQUEST_NULL 0
#  endif
#  ifndef MPI_MIN
#    define MPI_MIN 0
#  endif
#  ifndef MPI_MAX
#    define MPI_MAX 0
#  endif
#  ifndef MPI_SUM
#    define MPI_SUM 0
#  endif
#  ifndef MPI_LOR
#    define MPI_LOR 0
#  endif
#endif



/**
 * Helper macro to remove const from the pointer arguments to some MPI_*
 * functions.
 *
 * This is needed as the input arguments of functions like MPI_Allgather() are
 * not marked as const in OpenMPI 1.6.5. If using MPI 3 or newer, this macro
 * is a NOOP, while we do the following otherwise:
 *
 * 1. remove * from type of @p expr
 * 2. remove const from resulting type
 * 3. add * to resulting type
 * 4. const_cast the given expression @p expr to this new type.
 */
#ifdef DEAL_II_WITH_MPI
#  if DEAL_II_MPI_VERSION_GTE(3, 0)

#    define DEAL_II_MPI_CONST_CAST(expr) (expr)

#  else

#    include <type_traits>

#    define DEAL_II_MPI_CONST_CAST(expr)     \
      const_cast<typename std::remove_const< \
        typename std::remove_pointer<decltype(expr)>::type>::type *>(expr)

#  endif
#endif



DEAL_II_NAMESPACE_OPEN


// Forward type declarations to allow MPI sums over tensorial types
#ifndef DOXYGEN
template <int rank, int dim, typename Number>
class Tensor;
template <int rank, int dim, typename Number>
class SymmetricTensor;
template <typename Number>
class SparseMatrix;
class IndexSet;
#endif

namespace Utilities
{
  /**
   * Given the total number of elements @p total_size, create an evenly
   * distributed 1:1 partitioning of the elements for across @p n_partitions.
   * The local sizes will be equal to the @p total_size divided by the number
   * of partitions plus the remainder being divided amongst the first
   * processes. Each process will store a contiguous subset of indices, and the
   * index set on process p+1 starts at the index one larger than the last one
   * stored on process p.
   * For example, a @p total_size of 11 with 3 processes will result
   * in the IndexSets { [0,4), [4,8), [8,11)] }, and this function will
   * return the @p my_partition_id 's IndexSet.
   */
  IndexSet
  create_evenly_distributed_partitioning(const unsigned int my_partition_id,
                                         const unsigned int n_partitions,
                                         const IndexSet::size_type total_size);

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
     * @ref GlossMPICommunicator "communicator"
     * object. If this is a sequential job (i.e., the program
     * is not using MPI at all, or is using MPI but has been started with
     * only one MPI process), then the communicator necessarily involves
     * only one process and the function returns 1.
     */
    unsigned int
    n_mpi_processes(const MPI_Comm &mpi_communicator);

    /**
     * Return the
     * @ref GlossMPIRank "rank of the present MPI process"
     * in the space of processes described by the given
     * @ref GlossMPICommunicator "communicator".
     * This will be a unique value for each process between zero and (less
     * than) the number of all processes (given by get_n_mpi_processes()).
     */
    unsigned int
    this_mpi_process(const MPI_Comm &mpi_communicator);

    /**
     * Return a vector of the ranks (within @p comm_large) of a subset of
     * processes specified by @p comm_small.
     */
    const std::vector<unsigned int>
    mpi_processes_within_communicator(const MPI_Comm &comm_large,
                                      const MPI_Comm &comm_small);

    /**
     * Consider an unstructured communication pattern where every process in
     * an MPI universe wants to send some data to a subset of the other
     * processors. To do that, the other processors need to know who to expect
     * messages from. This function computes this information.
     *
     * @param mpi_comm A
     * @ref GlossMPICommunicator "communicator"
     * that describes the processors that are going to communicate with each
     * other.
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
    compute_point_to_point_communication_pattern(
      const MPI_Comm &                 mpi_comm,
      const std::vector<unsigned int> &destinations);

    /**
     * Simplified (for efficiency) version of the
     * compute_point_to_point_communication_pattern()
     * which only computes the number of processes in an MPI universe to expect
     * communication from.
     *
     * @param mpi_comm A
     * @ref GlossMPICommunicator "communicator"
     * that describes the processors that are going to communicate with each
     * other.
     *
     * @param destinations The list of processors the current process wants to
     * send information to. This list need not be sorted in any way. If it
     * contains duplicate entries that means that multiple messages are
     * intended for a given destination.
     *
     * @return A number of processors that want to send something to the current
     * processor.
     */
    unsigned int
    compute_n_point_to_point_communications(
      const MPI_Comm &                 mpi_comm,
      const std::vector<unsigned int> &destinations);

    /**
     * Given a
     * @ref GlossMPICommunicator "communicator",
     * generate a new communicator that contains the same set of processors
     * but that has a different, unique identifier.
     *
     * This functionality can be used to ensure that different objects, such
     * as distributed matrices, each have unique communicators over which they
     * can interact without interfering with each other.
     *
     * When no longer needed, the communicator created here needs to be
     * destroyed using free_communicator().
     *
     * This function is equivalent to calling
     * <code>MPI_Comm_dup(mpi_communicator, &return_value);</code>.
     */
    MPI_Comm
    duplicate_communicator(const MPI_Comm &mpi_communicator);

    /**
     * Free the given
     * @ref GlossMPICommunicator "communicator"
     * @p mpi_communicator that was duplicated using duplicate_communicator().
     *
     * The argument is passed by reference and will be invalidated and set to
     * the MPI null handle. This function is equivalent to calling
     * <code>MPI_Comm_free(&mpi_communicator);</code>.
     */
    void
    free_communicator(MPI_Comm &mpi_communicator);

    /**
     * Helper class to automatically duplicate and free an MPI
     * @ref GlossMPICommunicator "communicator".
     *
     * This class duplicates the communicator given in the constructor
     * using duplicate_communicator() and frees it automatically when
     * this object gets destroyed by calling free_communicator(). You
     * can access the wrapped communicator using operator*.
     *
     * This class exists to easily allow duplicating communicators without
     * having to worry when and how to free it after usage.
     */
    class DuplicatedCommunicator
    {
    public:
      /**
       * Create a duplicate of the given @p communicator.
       */
      explicit DuplicatedCommunicator(const MPI_Comm &communicator)
        : comm(duplicate_communicator(communicator))
      {}

      /**
       * Do not allow making copies.
       */
      DuplicatedCommunicator(const DuplicatedCommunicator &) = delete;

      /**
       * The destructor will free the communicator automatically.
       */
      ~DuplicatedCommunicator()
      {
        free_communicator(comm);
      }

      /**
       * Access the stored communicator.
       */
      const MPI_Comm &operator*() const
      {
        return comm;
      }


      /**
       * Do not allow assignment of this class.
       */
      DuplicatedCommunicator &
      operator=(const DuplicatedCommunicator &) = delete;

    private:
      /**
       * The communicator of course.
       */
      MPI_Comm comm;
    };

    /**
     * This class represents a mutex to guard a critical section for a set of
     * processors in a parallel computation using MPI.
     *
     * The lock() commands waits until all MPI ranks in the communicator have
     * released a previous lock using unlock().
     *
     * A typical usage involves guarding a critical section using a lock guard:
     * @code
     * {
     *   static CollectiveMutex      mutex;
     *   CollectiveMutex::ScopedLock lock(mutex, comm);
     *   // [ critical code to be guarded]
     * }
     * @endcode
     *
     * Here, the critical code will finish on all processors before the mutex
     * can be acquired again (for example by a second execution of the block
     * above. The critical code block typically involves MPI communication that
     * would yield incorrect results without the lock. For example, if the code
     * contains nonblocking receives with MPI_ANY_SOURCE, packets can be
     * confused between iterations.
     *
     * Note that the mutex needs to be the same instance between calls to the
     * same critical region. While not required, this can be achieved by making
     * the instance static (like in the example above). The variable can also be
     * a global variable, or a member variable of the object to which the
     * executing function belongs.
     */
    class CollectiveMutex
    {
    public:
      /**
       * This helper class provides a scoped lock for the CollectiveMutex.
       *
       * See the class documentation of CollectiveMutex for details.
       */
      class ScopedLock
      {
      public:
        /**
         * Constructor. Blocks until it can acquire the lock.
         */
        explicit ScopedLock(CollectiveMutex &mutex, const MPI_Comm &comm)
          : mutex(mutex)
          , comm(comm)
        {
          mutex.lock(comm);
        }

        /**
         * Destructor. Releases the lock.
         */
        ~ScopedLock()
        {
          mutex.unlock(comm);
        }

      private:
        /**
         * A reference to the mutex.
         */
        CollectiveMutex &mutex;
        /**
         * The communicator.
         */
        const MPI_Comm comm;
      };

      /**
       * Constructor of this class.
       */
      explicit CollectiveMutex();

      /**
       * Destroy the mutex. Assumes the lock is not currently held.
       */
      ~CollectiveMutex();

      /**
       * Acquire the mutex and, if necessary, wait until we can do so.
       *
       * This is a collective call that needs to be executed by all processors
       * in the communicator.
       */
      void
      lock(const MPI_Comm &comm);

      /**
       * Release the lock.
       *
       * This is a collective call that needs to be executed by all processors
       * in the communicator.
       */
      void
      unlock(const MPI_Comm &comm);

    private:
      /**
       * Keep track if we have this lock right now.
       */
      bool locked;

      /**
       * The request to keep track of the non-blocking barrier.
       */
      MPI_Request request;
    };



    /**
     * If @p comm is an intracommunicator, this function returns a new
     * communicator @p newcomm with communication group defined by the
     * @p group argument. The function is only collective over the group of
     * processes that actually want to create the communicator, i.e., that
     * are named in the @p group argument. If multiple threads at a given
     * process perform concurrent create_group() operations, the user must
     * distinguish these operations by providing different @p tag or @p comm
     * arguments.
     *
     * This function was introduced in the MPI-3.0 standard. If available,
     * the corresponding function in the provided MPI implementation is used.
     * Otherwise, the implementation follows the one described in the
     * following publication:
     * @code{.bib}
     * @inproceedings{dinan2011noncollective,
     *   title        = {Noncollective communicator creation in MPI},
     *   author       = {Dinan, James and Krishnamoorthy, Sriram and Balaji,
     *                   Pavan and Hammond, Jeff R and Krishnan, Manojkumar and
     *                   Tipparaju, Vinod and Vishnu, Abhinav},
     *   booktitle    = {European MPI Users' Group Meeting},
     *   pages        = {282--291},
     *   year         = {2011},
     *   organization = {Springer}
     * }
     * @endcode
     */
#ifdef DEAL_II_WITH_MPI
    int
    create_group(const MPI_Comm & comm,
                 const MPI_Group &group,
                 const int        tag,
                 MPI_Comm *       new_comm);
#endif

    /**
     * Given the number of locally owned elements @p locally_owned_size,
     * create a 1:1 partitioning of the of elements across the MPI
     * communicator @p comm. The total size of elements is the sum of
     * @p locally_owned_size across the MPI communicator.  Each process will
     * store contiguous subset of indices, and the index set on process p+1
     * starts at the index one larger than the last one stored on process p.
     */
    std::vector<IndexSet>
    create_ascending_partitioning(const MPI_Comm &          comm,
                                  const IndexSet::size_type locally_owned_size);

    /**
     * Given the total number of elements @p total_size, create an evenly
     * distributed 1:1 partitioning of the elements across the
     * MPI communicator @p comm.
     * Uses @p comm to determine number of partitions and processor ID to call the
     * @p create_evenly_distributed_partitioning() function above.
     */
    IndexSet
    create_evenly_distributed_partitioning(
      const MPI_Comm &          comm,
      const IndexSet::size_type total_size);

#ifdef DEAL_II_WITH_MPI
    /**
     * Calculate mean and standard deviation across the MPI communicator @p comm
     * for values provided as a range `[begin,end)`.
     * The mean is computed as $\bar x=\frac 1N \sum x_k$ where the $x_k$ are
     * the elements pointed to by the `begin` and `end` iterators on all
     * processors (i.e., each processor's `[begin,end)` range points to a subset
     * of the overall number of elements). The standard deviation is calculated
     * as $\sigma=\sqrt{\frac {1}{N-1} \sum |x_k -\bar x|^2}$, which is known as
     * unbiased sample variance.
     *
     * @tparam Number specifies the type to store the mean value.
     * The standard deviation is stored as the corresponding real type.
     * This allows, for example, to calculate statistics from integer input
     * values.
     */
    template <class Iterator, typename Number = long double>
    std::pair<Number, typename numbers::NumberTraits<Number>::real_type>
    mean_and_standard_deviation(const Iterator  begin,
                                const Iterator  end,
                                const MPI_Comm &comm);
#endif

    /**
     * Return the sum over all processors of the value @p t. This function is
     * collective over all processors given in the
     * @ref GlossMPICommunicator "communicator".
     * If deal.II is not configured for use of MPI, this function simply
     * returns the value of @p t. This function corresponds to the
     * <code>MPI_Allreduce</code> function, i.e. all processors receive the
     * result of this operation.
     *
     * @note Sometimes, not all processors need a result and in that case one
     * would call the <code>MPI_Reduce</code> function instead of the
     * <code>MPI_Allreduce</code> function. The latter is at most twice as
     * expensive, so if you are concerned about performance, it may be
     * worthwhile investigating whether your algorithm indeed needs the result
     * everywhere.
     *
     * @note This function is only implemented for certain template arguments
     * <code>T</code>, namely <code>float, double, int, unsigned int</code>.
     */
    template <typename T>
    T
    sum(const T &t, const MPI_Comm &mpi_communicator);

    /**
     * Like the previous function, but take the sums over the elements of an
     * array of type T. In other words, the i-th element of the results
     * array is the sum over the i-th entries of the input arrays from each
     * processor. T and U must decay to the same type, e.g. they just differ by
     * one of them having a const type qualifier and the other not.
     *
     * Input and output arrays may be the same.
     */
    template <typename T, typename U>
    void
    sum(const T &values, const MPI_Comm &mpi_communicator, U &sums);

    /**
     * Like the previous function, but take the sums over the elements of an
     * array as specified by the ArrayView arguments.
     * In other words, the i-th element of the results
     * array is the sum over the i-th entries of the input arrays from each
     * processor.
     *
     * Input and output arrays may be the same.
     */
    template <typename T>
    void
    sum(const ArrayView<const T> &values,
        const MPI_Comm &          mpi_communicator,
        const ArrayView<T> &      sums);

    /**
     * Perform an MPI sum of the entries of a symmetric tensor.
     *
     * @relatesalso SymmetricTensor
     */
    template <int rank, int dim, typename Number>
    SymmetricTensor<rank, dim, Number>
    sum(const SymmetricTensor<rank, dim, Number> &local,
        const MPI_Comm &                          mpi_communicator);

    /**
     * Perform an MPI sum of the entries of a tensor.
     *
     * @relatesalso Tensor
     */
    template <int rank, int dim, typename Number>
    Tensor<rank, dim, Number>
    sum(const Tensor<rank, dim, Number> &local,
        const MPI_Comm &                 mpi_communicator);

    /**
     * Perform an MPI sum of the entries of a SparseMatrix.
     *
     * @note @p local and @p global should have the same sparsity
     * pattern and it should be the same for all MPI processes.
     *
     * @relatesalso SparseMatrix
     */
    template <typename Number>
    void
    sum(const SparseMatrix<Number> &local,
        const MPI_Comm &            mpi_communicator,
        SparseMatrix<Number> &      global);

    /**
     * Return the maximum over all processors of the value @p t. This function
     * is collective over all processors given in the
     * @ref GlossMPICommunicator "communicator".
     * If deal.II is not configured for use of MPI, this function simply
     * returns the value of @p t. This function corresponds to the
     * <code>MPI_Allreduce</code> function, i.e. all processors receive the
     * result of this operation.
     *
     * @note Sometimes, not all processors need a result and in that case one
     * would call the <code>MPI_Reduce</code> function instead of the
     * <code>MPI_Allreduce</code> function. The latter is at most twice as
     * expensive, so if you are concerned about performance, it may be
     * worthwhile investigating whether your algorithm indeed needs the result
     * everywhere.
     *
     * @note This function is only implemented for certain template arguments
     * <code>T</code>, namely <code>float, double, int, unsigned int</code>.
     */
    template <typename T>
    T
    max(const T &t, const MPI_Comm &mpi_communicator);

    /**
     * Like the previous function, but take the maximum over the elements of an
     * array of type T. In other words, the i-th element of the results array is
     * the maximum over the i-th entries of the input arrays from each
     * processor. T and U must decay to the same type, e.g. they just differ by
     * one of them having a const type qualifier and the other not.
     *
     * Input and output vectors may be the same.
     */
    template <typename T, typename U>
    void
    max(const T &values, const MPI_Comm &mpi_communicator, U &maxima);

    /**
     * Like the previous function, but take the maximum over the elements of an
     * array as specified by the ArrayView arguments.
     * In other words, the i-th element of the results
     * array is the maximum over the i-th entries of the input arrays from each
     * processor.
     *
     * Input and output arrays may be the same.
     */
    template <typename T>
    void
    max(const ArrayView<const T> &values,
        const MPI_Comm &          mpi_communicator,
        const ArrayView<T> &      maxima);

    /**
     * Return the minimum over all processors of the value @p t. This function
     * is collective over all processors given in the
     * @ref GlossMPICommunicator "communicator".
     * If deal.II is not configured for use of MPI, this function simply
     * returns the value of @p t. This function corresponds to the
     * <code>MPI_Allreduce</code> function, i.e. all processors receive the
     * result of this operation.
     *
     * @note Sometimes, not all processors need a result and in that case one
     * would call the <code>MPI_Reduce</code> function instead of the
     * <code>MPI_Allreduce</code> function. The latter is at most twice as
     * expensive, so if you are concerned about performance, it may be
     * worthwhile investigating whether your algorithm indeed needs the result
     * everywhere.
     *
     * @note This function is only implemented for certain template arguments
     * <code>T</code>, namely <code>float, double, int, unsigned int</code>.
     */
    template <typename T>
    T
    min(const T &t, const MPI_Comm &mpi_communicator);

    /**
     * Like the previous function, but take the minima over the elements of an
     * array of type T. In other words, the i-th element of the results
     * array is the minimum of the i-th entries of the input arrays from each
     * processor. T and U must decay to the same type, e.g. they just differ by
     * one of them having a const type qualifier and the other not.
     *
     * Input and output arrays may be the same.
     */
    template <typename T, typename U>
    void
    min(const T &values, const MPI_Comm &mpi_communicator, U &minima);

    /**
     * Like the previous function, but take the minimum over the elements of an
     * array as specified by the ArrayView arguments.
     * In other words, the i-th element of the results
     * array is the minimum over the i-th entries of the input arrays from each
     * processor.
     *
     * Input and output arrays may be the same.
     */
    template <typename T>
    void
    min(const ArrayView<const T> &values,
        const MPI_Comm &          mpi_communicator,
        const ArrayView<T> &      minima);

    /**
     * Performs a <i>logical or</i> operation over all processors of the value
     * @p t. The <i>logical or</i> operator `||` returns the boolean value
     * `true` if either or all operands are `true` and returns `false`
     * otherwise. If the provided value @p t corresponds to `0` in its
     * associated data type `T`, it will be interpreted as `false`, and `true`
     * otherwise. Data type `T` must be of type `integral`, i.e., `bool`,
     * `char`, `short`, `int`, `long`, or any of their variations.
     *
     * This function is collective over all processors given in the
     * @ref GlossMPICommunicator "communicator".
     * If deal.II is not configured for use of MPI, this function simply
     * returns the value of @p value. This function corresponds to the
     * <code>MPI_Allreduce</code> function, i.e., all processors receive the
     * result of this operation.
     *
     * @note Sometimes, not all processors need a result and in that case one
     * would call the <code>MPI_Reduce</code> function instead of the
     * <code>MPI_Allreduce</code> function. The latter is at most twice as
     * expensive, so if you are concerned about performance, it may be
     * worthwhile investigating whether your algorithm indeed needs the result
     * everywhere.
     */
    template <typename T>
    T
    logical_or(const T &t, const MPI_Comm &mpi_communicator);

    /**
     * Like the previous function, but performs the <i>logical or</i> operation
     * on each element of an array. In other words, the i-th element of the
     * results array is the result of the <i>logical or</i> operation applied on
     * the i-th entries of the input arrays from each processor. T and U must
     * decay to the same type, e.g., they just differ by one of them having a
     * const type qualifier and the other not.
     *
     * Input and output arrays may be the same.
     *
     * @note Depending on your standard library, this function may not work with
     *   specializations of `std::vector` for the data type `bool`. In that
     *   case, use a different container or data type.
     */
    template <typename T, typename U>
    void
    logical_or(const T &values, const MPI_Comm &mpi_communicator, U &results);

    /**
     * Like the previous function, but performs the <i>logical or</i> operation
     * on each element of an array as specified by the ArrayView arguments.
     * In other words, the i-th element of the results array is the result of
     * the <i>logical or</i> operation applied on the i-th entries of the input
     * arrays from each processor.
     *
     * Input and output arrays may be the same.
     */
    template <typename T>
    void
    logical_or(const ArrayView<const T> &values,
               const MPI_Comm &          mpi_communicator,
               const ArrayView<T> &      results);

    /**
     * A data structure to store the result of the min_max_avg() function.
     * The structure stores the minimum, maximum, and average of one
     * value contributed by each processor that participates in an
     * @ref GlossMPICommunicator "MPI communicator".
     * The structure also stores
     * the indices (or, more precisely, the
     * @ref GlossMPIRank "MPI rank")
     * of the processors that hold the minimum and maximum values,
     * as well as the sum over all values.
     *
     * @note This structure has no constructors because MPI requires it
     *   to be a POD type.
     */
    struct MinMaxAvg
    {
      /**
       * The sum over all values contributed by the processors that
       * participate in the call to min_max_avg().
       */
      double sum;

      /**
       * The minimum value over all values contributed by the processors that
       * participate in the call to min_max_avg().
       */
      double min;

      /**
       * The maximum value over all values contributed by the processors that
       * participate in the call to min_max_avg().
       */
      double max;

      /**
       * One of the ranks (i.e.,
       * @ref GlossMPIRank "MPI rank"
       * within an
       * @ref GlossMPICommunicator "MPI communicator")
       * of the
       * processors that hold the minimal value.
       */
      unsigned int min_index;

      /**
       * One of the ranks (i.e.,
       * @ref GlossMPIRank "MPI rank"
       * within an
       * @ref GlossMPICommunicator "MPI communicator")
       * of the
       * processors that hold the maximal value.
       */
      unsigned int max_index;

      /**
       * The average of the values contributed by the processors that
       * participate in the call to min_max_avg().
       */
      double avg;
    };

    /**
     * Return sum, average, minimum, maximum, processor id of minimum and
     * maximum as a collective operation of on the given MPI
     * @ref GlossMPICommunicator "communicator"
     * @p mpi_communicator. Each processor's value is given in @p my_value and
     * the result will be returned. The result is available on all machines.
     *
     * @note Sometimes, not all processors need a result and in that case one
     * would call the <code>MPI_Reduce</code> function instead of the
     * <code>MPI_Allreduce</code> function. The latter is at most twice as
     * expensive, so if you are concerned about performance, it may be
     * worthwhile investigating whether your algorithm indeed needs the result
     * everywhere.
     */
    MinMaxAvg
    min_max_avg(const double my_value, const MPI_Comm &mpi_communicator);

    /**
     * Same as above but returning the sum, average, minimum, maximum,
     * process id of minimum and maximum as a collective operation on the
     * given MPI
     * @ref GlossMPICommunicator "communicator"
     * @p mpi_communicator for each entry of the vector.
     *
     * @note This function performs a single reduction sweep.
     *
     * @pre Size of the input vector has to be the same on all processes.
     */
    std::vector<MinMaxAvg>
    min_max_avg(const std::vector<double> &my_value,
                const MPI_Comm &           mpi_communicator);


    /**
     * Same as above but returning the sum, average, minimum, maximum,
     * process id of minimum and maximum as a collective operation on the
     * given MPI
     * @ref GlossMPICommunicator "communicator"
     * @p mpi_communicator for each entry of the ArrayView.
     *
     * @note This function performs a single reduction sweep.
     *
     * @pre Size of the input ArrayView has to be the same on all processes
     *   and the input and output ArrayVew have to have the same size.
     */
    void
    min_max_avg(const ArrayView<const double> &my_values,
                const ArrayView<MinMaxAvg> &   result,
                const MPI_Comm &               mpi_communicator);


    /**
     * A class that is used to initialize the MPI system at the beginning of a
     * program and to shut it down again at the end. It also allows you to
     * control the number of threads used within each MPI process.
     *
     * If deal.II is configured with PETSc, PETSc will be initialized
     * via `PetscInitialize` in the beginning (constructor of this
     * class) and de-initialized via `PetscFinalize` at the end (i.e.,
     * in the destructor of this class). The same is true for SLEPc.
     *
     * If deal.II is configured with p4est, that library will also be
     * initialized in the beginning, and de-initialized at the end
     * (by calling sc_init(), p4est_init(), and sc_finalize()).
     *
     * If a program uses MPI one would typically just create an object
     * of this type at the beginning of <code>main()</code>. The
     * constructor of this class then runs <code>MPI_Init()</code>
     * with the given arguments and also initializes the other
     * libraries mentioned above. At the end of the program, the
     * compiler will invoke the destructor of this object which in
     * turns calls <code>MPI_Finalize</code> to shut down the MPI
     * system.
     *
     * This class is used in step-17, step-18, step-40, step-32, and
     * several others.
     *
     * @note This class performs initialization of the MPI subsystem
     * as well as the dependent libraries listed above through the
     * `MPI_COMM_WORLD` communicator. This means that you will have to
     * create an MPI_InitFinalize object on <i>all</i> MPI processes,
     * whether or not you intend to use deal.II on a given
     * processor. In most use cases, one will of course want to work
     * on all MPI processes using essentially the same program, and so
     * this is not an issue. But if you plan to run deal.II-based work
     * on only a subset of MPI processes, using an @ ref
     * GlossMPICommunicator "MPI communicator" that is a subset of
     * `MPI_COMM_WORLD` (for example, in client-server settings where
     * only a subset of processes is responsible for the finite
     * element communications and the remaining processes do other
     * things), then you still need to create this object here on all
     * MPI processes at the beginning of the program because it uses
     * `MPI_COMM_WORLD` during initialization.
     */
    class MPI_InitFinalize
    {
    public:
      /**
       * Initialize MPI (and, if deal.II was configured to use it, PETSc) and
       * set the number of threads used by deal.II (via the underlying
       * Threading Building Blocks library) to the given parameter.
       *
       * @param[in,out] argc A reference to the 'argc' argument passed to
       * main. This argument is used to initialize MPI (and, possibly, PETSc)
       * as they read arguments from the command line.
       * @param[in,out] argv A reference to the 'argv' argument passed to
       * main.
       * @param[in] max_num_threads The maximal number of threads this MPI
       * process should utilize. If this argument is set to
       * numbers::invalid_unsigned_int (the default value), then the number of
       * threads is determined automatically in the following way: the number
       * of threads to run on this MPI process is set in such a way that all
       * of the cores in your node are spoken for. In other words, if you have
       * started one MPI process per node, setting this argument is equivalent
       * to setting it to the number of cores present in the node this MPI
       * process runs on. If you have started as many MPI processes per node
       * as there are cores on each node, then this is equivalent to passing 1
       * as the argument. On the other hand, if, for example, you start 4 MPI
       * processes on each 16-core node, then this option will start 4 worker
       * threads for each node. If you start 3 processes on an 8 core node,
       * then they will start 3, 3 and 2 threads, respectively.
       *
       * @note This function calls MultithreadInfo::set_thread_limit() with
       * either @p max_num_threads or, following the discussion above, a
       * number of threads equal to the number of cores allocated to this MPI
       * process. However, MultithreadInfo::set_thread_limit() in turn also
       * evaluates the environment variable DEAL_II_NUM_THREADS. Finally, the
       * worker threads can only be created on cores to which the current MPI
       * process has access to; some MPI implementations limit the number of
       * cores each process may access to one or a subset of cores in order to
       * ensure better cache behavior. Consequently, the number of threads
       * that will really be created will be the minimum of the argument
       * passed here, the environment variable (if set), and the number of
       * cores accessible to the thread.
       *
       * @note MultithreadInfo::set_thread_limit() can only work if it is
       * called before any threads are created. The safest place for a call to
       * it is therefore at the beginning of <code>main()</code>.
       * Consequently, this extends to the current class: the best place to
       * create an object of this type is also at or close to the top of
       * <code>main()</code>.
       */
      MPI_InitFinalize(
        int &              argc,
        char **&           argv,
        const unsigned int max_num_threads = numbers::invalid_unsigned_int);

      /**
       * Destructor. Calls <tt>MPI_Finalize()</tt> in case this class owns the
       * MPI process.
       */
      ~MPI_InitFinalize();

      /**
       * Register a reference to an MPI_Request
       * on which we need to call `MPI_Wait` before calling `MPI_Finalize`.
       *
       * The object @p request needs to exist when MPI_Finalize is called, which means the
       * request is typically statically allocated. Otherwise, you need to call
       * unregister_request() before the request goes out of scope. Note that it
       * is acceptable for a request to be already waited on (and consequently
       * reset to MPI_REQUEST_NULL).
       *
       * It is acceptable to call this function more than once with the same
       * instance (as it is done in the example below).
       *
       * Typically, this function is used by CollectiveMutex and not directly,
       * but it can also be used directly like this:
       * @code
       * void my_fancy_communication()
       * {
       *   static MPI_Request request = MPI_REQUEST_NULL;
       *   MPI_InitFinalize::register_request(request);
       *   MPI_Wait(&request, MPI_STATUS_IGNORE);
       *   // [some algorithm that is not safe to be executed twice in a row.]
       *   MPI_IBarrier(comm, &request);
       * }
       * @endcode
       */
      static void
      register_request(MPI_Request &request);

      /**
       * Unregister a request previously added using register_request().
       */
      static void
      unregister_request(MPI_Request &request);

      /**
       * A structure that has boost::signal objects to register a call back
       * to run after MPI init or finalize.
       *
       * For documentation on signals, see
       * http://www.boost.org/doc/libs/release/libs/signals2 .
       */
      struct Signals
      {
        /**
         * A signal that is triggered immediately after we have
         * initialized the MPI context with <code>MPI_Init()</code>.
         */
        boost::signals2::signal<void()> at_mpi_init;

        /**
         * A signal that is triggered just before we close the MPI context
         * with <code>MPI_Finalize()</code>. It can be used to deallocate
         * statically allocated MPI resources that need to be deallocated
         * before <code>MPI_Finalize()</code> is called.
         */
        boost::signals2::signal<void()> at_mpi_finalize;
      };

      static Signals signals;

    private:
      /**
       * Requests to MPI_Wait before finalizing
       */
      static std::set<MPI_Request *> requests;
    };

    /**
     * Return whether (i) deal.II has been compiled to support MPI (for
     * example by compiling with <code>CXX=mpiCC</code>) and if so whether
     * (ii) <code>MPI_Init()</code> has been called (for example using the
     * Utilities::MPI::MPI_InitFinalize class). In other words, the result
     * indicates whether the current job is running under MPI.
     *
     * @note The function does not take into account whether an MPI job
     * actually runs on more than one processor or is, in fact, a single-node
     * job that happens to run under MPI.
     */
    bool
    job_supports_mpi();

    /**
     * Initiate a some-to-some communication, and exchange arbitrary objects
     * (the class T should be serializable using boost::serialize) between
     * processors.
     *
     * @param[in] comm MPI communicator.
     *
     * @param[in] objects_to_send A map from the rank (unsigned int) of the
     *  process meant to receive the data and the object to send (the type `T`
     *  must be serializable for this function to work properly). If this
     *  map contains an entry with a key equal to the rank of the current
     *  process (i.e., an instruction to a process to send data to itself),
     *  then this data item is simply copied to the returned object.
     *
     * @return A map from the rank (unsigned int) of the process
     *  which sent the data and object received.
     */
    template <typename T>
    std::map<unsigned int, T>
    some_to_some(const MPI_Comm &                 comm,
                 const std::map<unsigned int, T> &objects_to_send);

    /**
     * A generalization of the classic MPI_Allgather function, that accepts
     * arbitrary data types T, as long as boost::serialize accepts T as an
     * argument.
     *
     * @param[in] comm MPI communicator.
     * @param[in] object_to_send An object to send to all other processes
     *
     * @return A vector of objects, with size equal to the number of
     *  processes in the MPI communicator. Each entry contains the object
     *  received from the processor with the corresponding rank within the
     *  communicator.
     */
    template <typename T>
    std::vector<T>
    all_gather(const MPI_Comm &comm, const T &object_to_send);

    /**
     * A generalization of the classic MPI_Gather function, that accepts
     * arbitrary data types T, as long as boost::serialize accepts T as an
     * argument.
     *
     * @param[in] comm MPI communicator.
     * @param[in] object_to_send an object to send to the root process
     * @param[in] root_process The process, which receives the objects from all
     * processes. By default the process with rank 0 is the root process.
     *
     * @return The @p root_process receives a vector of objects, with size equal to the number of
     *  processes in the MPI communicator. Each entry contains the object
     *  received from the processor with the corresponding rank within the
     *  communicator. All other processes receive an empty vector.
     */
    template <typename T>
    std::vector<T>
    gather(const MPI_Comm &   comm,
           const T &          object_to_send,
           const unsigned int root_process = 0);

    /**
     * Sends an object @p object_to_send from the process @p root_process
     * to all other processes.
     *
     * A generalization of the classic `MPI_Bcast` function that accepts
     * arbitrary data types `T`, as long as Utilities::pack() (which in turn
     * uses `boost::serialize`, see in Utilities::pack() for details) accepts
     * `T` as an argument.
     *
     * @param[in] comm MPI communicator.
     * @param[in] object_to_send An object to send to all processes.
     * @param[in] root_process The process that sends the object to all
     * processes. By default the process with rank 0 is the root process.
     *
     * @return On the root process, return a copy of @p object_to_send.
     *   On every other process, return a copy of the object sent by
     *   the @p root_process.
     */
    template <typename T>
    T
    broadcast(const MPI_Comm &   comm,
              const T &          object_to_send,
              const unsigned int root_process = 0);

    /**
     * A function that combines values @p local_value from all processes
     * via a user-specified binary operation @p combiner and distributes the
     * result back to all processes. As such this function is similar to
     * MPI_Allreduce (if it were implemented by a global reduction followed
     * by a broadcast step) but due to the user-specified binary operation also
     * general object types, including ones that store variable amounts of data,
     * can be handled.
     */
    template <typename T>
    T
    all_reduce(const T &                                     local_value,
               const MPI_Comm &                              comm,
               const std::function<T(const T &, const T &)> &combiner);

    /**
     * Given a partitioned index set space, compute the owning MPI process rank
     * of each element of a second index set according to the partitioned index
     * set. A natural usage of this function is to compute for each ghosted
     * degree of freedom the MPI rank of the process owning that index.
     *
     * One might think: "But we know which rank a ghost DoF belongs to based on
     * the subdomain id of the cell it is on". But this heuristic fails for DoFs
     * on interfaces between ghost cells with different subdomain_ids, or
     * between a ghost cell and an artificial cell. Furthermore, this function
     * enables a completely abstract exchange of information without the help of
     * the mesh in terms of neighbors.
     *
     * The first argument passed to this function, @p owned_indices, must
     * uniquely partition an index space between all processes.
     * Otherwise, there are no limitations on this argument: In particular,
     * there is no need in partitioning
     * the index space into contiguous subsets. Furthermore, there are no
     * limitations
     * on the second index set @p indices_to_look_up as long as the size matches
     * the first one. It can be chosen arbitrarily and independently on each
     * process. In the case that the second index set also contains locally
     * owned indices, these indices will be treated correctly and the rank of
     * this process is returned for those entries.
     *
     * @note This is a collective operation: all processes within the given
     * communicator have to call this function. Since this function does not
     * use MPI_Alltoall or MPI_Allgather, but instead uses non-blocking
     * point-to-point communication instead, and only a single non-blocking
     * barrier, it reduces the memory consumption significantly. This function
     * is suited for large-scale simulations with >100k MPI ranks.
     *
     * @param[in] owned_indices Index set with indices locally owned by this
     *            process.
     * @param[in] indices_to_look_up Index set containing indices of which the
     *            user is interested the rank of the owning process.
     * @param[in] comm MPI communicator.
     *
     * @return List containing the MPI process rank for each entry in the index
     *         set @p indices_to_look_up. The order coincides with the order
     *         within the ElementIterator.
     */
    std::vector<unsigned int>
    compute_index_owner(const IndexSet &owned_indices,
                        const IndexSet &indices_to_look_up,
                        const MPI_Comm &comm);

    /**
     * Compute the union of the input vectors @p vec of all processes in the
     *   MPI communicator @p comm.
     *
     * @note This is a collective operation. The result will available on all
     *   processes.
     */
    template <typename T>
    std::vector<T>
    compute_set_union(const std::vector<T> &vec, const MPI_Comm &comm);

    /**
     * The same as above but for std::set.
     */
    template <typename T>
    std::set<T>
    compute_set_union(const std::set<T> &set, const MPI_Comm &comm);

#ifndef DOXYGEN
    // declaration for an internal function that lives in mpi.templates.h
    namespace internal
    {
      template <typename T>
      void
      all_reduce(const MPI_Op &            mpi_op,
                 const ArrayView<const T> &values,
                 const MPI_Comm &          mpi_communicator,
                 const ArrayView<T> &      output);
    }

    // Since these depend on N they must live in the header file
    template <typename T, unsigned int N>
    void
    sum(const T (&values)[N], const MPI_Comm &mpi_communicator, T (&sums)[N])
    {
      internal::all_reduce(MPI_SUM,
                           ArrayView<const T>(values, N),
                           mpi_communicator,
                           ArrayView<T>(sums, N));
    }

    template <typename T, unsigned int N>
    void
    max(const T (&values)[N], const MPI_Comm &mpi_communicator, T (&maxima)[N])
    {
      internal::all_reduce(MPI_MAX,
                           ArrayView<const T>(values, N),
                           mpi_communicator,
                           ArrayView<T>(maxima, N));
    }

    template <typename T, unsigned int N>
    void
    min(const T (&values)[N], const MPI_Comm &mpi_communicator, T (&minima)[N])
    {
      internal::all_reduce(MPI_MIN,
                           ArrayView<const T>(values, N),
                           mpi_communicator,
                           ArrayView<T>(minima, N));
    }

    template <typename T, unsigned int N>
    void
    logical_or(const T (&values)[N],
               const MPI_Comm &mpi_communicator,
               T (&results)[N])
    {
      static_assert(std::is_integral<T>::value,
                    "The MPI_LOR operation only allows integral data types.");

      internal::all_reduce(MPI_LOR,
                           ArrayView<const T>(values, N),
                           mpi_communicator,
                           ArrayView<T>(results, N));
    }

    template <typename T>
    std::map<unsigned int, T>
    some_to_some(const MPI_Comm &                 comm,
                 const std::map<unsigned int, T> &objects_to_send)
    {
#  ifndef DEAL_II_WITH_MPI
      (void)comm;
      Assert(objects_to_send.size() < 2,
             ExcMessage("Cannot send to more than one processor."));
      Assert(objects_to_send.find(0) != objects_to_send.end() ||
               objects_to_send.size() == 0,
             ExcMessage("Can only send to myself or to nobody."));
      return objects_to_send;
#  else
      const auto my_proc = this_mpi_process(comm);

      std::map<unsigned int, T> received_objects;

      std::vector<unsigned int> send_to;
      send_to.reserve(objects_to_send.size());
      for (const auto &m : objects_to_send)
        if (m.first == my_proc)
          received_objects[my_proc] = m.second;
        else
          send_to.emplace_back(m.first);

      const unsigned int n_point_point_communications =
        Utilities::MPI::compute_n_point_to_point_communications(comm, send_to);

      // Protect the following communication:
      static CollectiveMutex      mutex;
      CollectiveMutex::ScopedLock lock(mutex, comm);

      // If we have something to send, or we expect something from other
      // processors, we need to visit one of the two scopes below. Otherwise,
      // no other action is required by this mpi process, and we can safely
      // return.
      if (send_to.size() == 0 && n_point_point_communications == 0)
        return received_objects;

      const int mpi_tag =
        internal::Tags::compute_point_to_point_communication_pattern;

      // Sending buffers
      std::vector<std::vector<char>> buffers_to_send(send_to.size());
      std::vector<MPI_Request>       buffer_send_requests(send_to.size());
      {
        unsigned int i = 0;
        for (const auto &rank_obj : objects_to_send)
          if (rank_obj.first != my_proc)
            {
              const auto &rank   = rank_obj.first;
              buffers_to_send[i] = Utilities::pack(rank_obj.second,
                                                   /*allow_compression=*/false);
              const int ierr     = MPI_Isend(buffers_to_send[i].data(),
                                         buffers_to_send[i].size(),
                                         MPI_CHAR,
                                         rank,
                                         mpi_tag,
                                         comm,
                                         &buffer_send_requests[i]);
              AssertThrowMPI(ierr);
              ++i;
            }
      }

      // Fill the output map
      {
        std::vector<char> buffer;
        // We do this on a first come/first served basis
        for (unsigned int i = 0; i < n_point_point_communications; ++i)
          {
            // Probe what's going on. Take data from the first available sender
            MPI_Status status;
            int        ierr = MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
            AssertThrowMPI(ierr);

            // Length of the message
            int len;
            ierr = MPI_Get_count(&status, MPI_CHAR, &len);
            AssertThrowMPI(ierr);
            buffer.resize(len);

            // Source rank
            const unsigned int rank = status.MPI_SOURCE;

            // Actually receive the message
            ierr = MPI_Recv(buffer.data(),
                            len,
                            MPI_CHAR,
                            status.MPI_SOURCE,
                            status.MPI_TAG,
                            comm,
                            MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);
            Assert(received_objects.find(rank) == received_objects.end(),
                   ExcInternalError(
                     "I should not receive again from this rank"));
            received_objects[rank] =
              Utilities::unpack<T>(buffer,
                                   /*allow_compression=*/false);
          }
      }

      // Wait to have sent all objects.
      const int ierr = MPI_Waitall(send_to.size(),
                                   buffer_send_requests.data(),
                                   MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);

      return received_objects;
#  endif // deal.II with MPI
    }

    template <typename T>
    std::vector<T>
    all_gather(const MPI_Comm &comm, const T &object)
    {
      if (job_supports_mpi() == false)
        return {object};

#  ifndef DEAL_II_WITH_MPI
      (void)comm;
      std::vector<T> v(1, object);
      return v;
#  else
      const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);

      std::vector<char> buffer = Utilities::pack(object);

      int n_local_data = buffer.size();

      // Vector to store the size of loc_data_array for every process
      std::vector<int> size_all_data(n_procs, 0);

      // Exchanging the size of each buffer
      MPI_Allgather(
        &n_local_data, 1, MPI_INT, size_all_data.data(), 1, MPI_INT, comm);

      // Now computing the displacement, relative to recvbuf,
      // at which to store the incoming buffer
      std::vector<int> rdispls(n_procs);
      rdispls[0] = 0;
      for (unsigned int i = 1; i < n_procs; ++i)
        rdispls[i] = rdispls[i - 1] + size_all_data[i - 1];

      // Step 3: exchange the buffer:
      std::vector<char> received_unrolled_buffer(rdispls.back() +
                                                 size_all_data.back());

      MPI_Allgatherv(buffer.data(),
                     n_local_data,
                     MPI_CHAR,
                     received_unrolled_buffer.data(),
                     size_all_data.data(),
                     rdispls.data(),
                     MPI_CHAR,
                     comm);

      std::vector<T> received_objects(n_procs);
      for (unsigned int i = 0; i < n_procs; ++i)
        {
          std::vector<char> local_buffer(received_unrolled_buffer.begin() +
                                           rdispls[i],
                                         received_unrolled_buffer.begin() +
                                           rdispls[i] + size_all_data[i]);
          received_objects[i] = Utilities::unpack<T>(local_buffer);
        }

      return received_objects;
#  endif
    }

    template <typename T>
    std::vector<T>
    gather(const MPI_Comm &   comm,
           const T &          object_to_send,
           const unsigned int root_process)
    {
#  ifndef DEAL_II_WITH_MPI
      (void)comm;
      (void)root_process;
      std::vector<T> v(1, object_to_send);
      return v;
#  else
      const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);
      const auto my_rank = dealii::Utilities::MPI::this_mpi_process(comm);

      AssertIndexRange(root_process, n_procs);

      std::vector<char> buffer       = Utilities::pack(object_to_send);
      int               n_local_data = buffer.size();

      // Vector to store the size of loc_data_array for every process
      // only the root process needs to allocate memory for that purpose
      std::vector<int> size_all_data;
      if (my_rank == root_process)
        size_all_data.resize(n_procs, 0);

      // Exchanging the size of each buffer
      int ierr = MPI_Gather(&n_local_data,
                            1,
                            MPI_INT,
                            size_all_data.data(),
                            1,
                            MPI_INT,
                            root_process,
                            comm);
      AssertThrowMPI(ierr);

      // Now computing the displacement, relative to recvbuf,
      // at which to store the incoming buffer; only for root
      std::vector<int> rdispls;
      if (my_rank == root_process)
        {
          rdispls.resize(n_procs, 0);
          for (unsigned int i = 1; i < n_procs; ++i)
            rdispls[i] = rdispls[i - 1] + size_all_data[i - 1];
        }
      // exchange the buffer:
      std::vector<char> received_unrolled_buffer;
      if (my_rank == root_process)
        received_unrolled_buffer.resize(rdispls.back() + size_all_data.back());

      ierr = MPI_Gatherv(buffer.data(),
                         n_local_data,
                         MPI_CHAR,
                         received_unrolled_buffer.data(),
                         size_all_data.data(),
                         rdispls.data(),
                         MPI_CHAR,
                         root_process,
                         comm);
      AssertThrowMPI(ierr);

      std::vector<T> received_objects;

      if (my_rank == root_process)
        {
          received_objects.resize(n_procs);

          for (unsigned int i = 0; i < n_procs; ++i)
            {
              const std::vector<char> local_buffer(
                received_unrolled_buffer.begin() + rdispls[i],
                received_unrolled_buffer.begin() + rdispls[i] +
                  size_all_data[i]);
              received_objects[i] = Utilities::unpack<T>(local_buffer);
            }
        }
      return received_objects;
#  endif
    }



    template <typename T>
    T
    broadcast(const MPI_Comm &   comm,
              const T &          object_to_send,
              const unsigned int root_process)
    {
#  ifndef DEAL_II_WITH_MPI
      (void)comm;
      (void)root_process;
      return object_to_send;
#  else
      const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);
      AssertIndexRange(root_process, n_procs);
      (void)n_procs;

      std::vector<char> buffer;
      unsigned int      buffer_size = numbers::invalid_unsigned_int;

      // On the root process, pack the data and determine what the
      // buffer size needs to be.
      if (this_mpi_process(comm) == root_process)
        {
          buffer      = Utilities::pack(object_to_send, false);
          buffer_size = buffer.size();
        }

      // Exchange the size of buffer
      int ierr = MPI_Bcast(&buffer_size, 1, MPI_UNSIGNED, root_process, comm);
      AssertThrowMPI(ierr);

      // If not on the root process, correctly size the buffer to
      // receive the data, then do exactly that.
      if (this_mpi_process(comm) != root_process)
        buffer.resize(buffer_size);

      ierr =
        MPI_Bcast(buffer.data(), buffer_size, MPI_CHAR, root_process, comm);
      AssertThrowMPI(ierr);

      if (Utilities::MPI::this_mpi_process(comm) == root_process)
        return object_to_send;
      else
        return Utilities::unpack<T>(buffer, false);
#  endif
    }


#  ifdef DEAL_II_WITH_MPI
    template <class Iterator, typename Number>
    std::pair<Number, typename numbers::NumberTraits<Number>::real_type>
    mean_and_standard_deviation(const Iterator  begin,
                                const Iterator  end,
                                const MPI_Comm &comm)
    {
      // below we do simple and straight-forward implementation. More elaborate
      // options are:
      // http://dx.doi.org/10.1145/2807591.2807644 section 3.1.2
      // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
      // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online
      using Std        = typename numbers::NumberTraits<Number>::real_type;
      const Number sum = std::accumulate(begin, end, Number(0.));

      const auto size = Utilities::MPI::sum(std::distance(begin, end), comm);
      Assert(size > 0, ExcDivideByZero());
      const Number mean =
        Utilities::MPI::sum(sum, comm) / static_cast<Std>(size);
      Std sq_sum = 0.;
      std::for_each(begin, end, [&mean, &sq_sum](const Number &v) {
        sq_sum += numbers::NumberTraits<Number>::abs_square(v - mean);
      });
      sq_sum = Utilities::MPI::sum(sq_sum, comm);
      return std::make_pair(mean,
                            std::sqrt(sq_sum / static_cast<Std>(size - 1)));
    }
#  endif

#endif
  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
