// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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

#ifndef dealii_mpi_h
#define dealii_mpi_h

#include <deal.II/base/config.h>
#include <deal.II/base/array_view.h>

#include <vector>
#include <map>

#if !defined(DEAL_II_WITH_MPI) && !defined(DEAL_II_WITH_PETSC)
// without MPI, we would still like to use
// some constructs with MPI data
// types. Therefore, create some dummies
typedef int MPI_Comm;
typedef int MPI_Datatype;
typedef int MPI_Op;
#  ifndef MPI_COMM_WORLD
#    define MPI_COMM_WORLD 0
#  endif
#  ifndef MPI_COMM_SELF
#    define MPI_COMM_SELF 0
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
#endif

DEAL_II_NAMESPACE_OPEN


//Forward type declarations to allow MPI sums over tensorial types
template <int rank, int dim, typename Number> class Tensor;
template <int rank, int dim, typename Number> class SymmetricTensor;
template <typename Number> class SparseMatrix;

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
     * @ref GlossMPICommunicator "communicator"
     * object. If this is a sequential job (i.e., the program
     * is not using MPI at all, or is using MPI but has been started with
     * only one MPI process), then the communicator necessarily involves
     * only one process and the function returns 1.
     */
    unsigned int n_mpi_processes (const MPI_Comm &mpi_communicator);

    /**
     * Return the
     * @ref GlossMPIRank "rank of the present MPI process"
     * in the space of processes described by the given
     * @ref GlossMPICommunicator "communicator".
     * This will be a unique value for each process between zero and (less
     * than) the number of all processes (given by get_n_mpi_processes()).
     */
    unsigned int this_mpi_process (const MPI_Comm &mpi_communicator);

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
    compute_point_to_point_communication_pattern (const MPI_Comm &mpi_comm,
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
     * destroyed using <code>MPI_Comm_free</code>.
     */
    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator);

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
    T sum (const T &t,
           const MPI_Comm &mpi_communicator);

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
    void sum (const T        &values,
              const MPI_Comm &mpi_communicator,
              U              &sums);

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
    void sum (const ArrayView<const T> &values,
              const MPI_Comm           &mpi_communicator,
              const ArrayView<T>       &sums);

    /**
     * Perform an MPI sum of the entries of a symmetric tensor.
     *
     * @relatesalso SymmetricTensor
     */
    template <int rank, int dim, typename Number>
    SymmetricTensor<rank,dim,Number>
    sum (const SymmetricTensor<rank,dim,Number> &local,
         const MPI_Comm &mpi_communicator);

    /**
     * Perform an MPI sum of the entries of a tensor.
     *
     * @relatesalso Tensor
     */
    template <int rank, int dim, typename Number>
    Tensor<rank,dim,Number>
    sum (const Tensor<rank,dim,Number> &local,
         const MPI_Comm &mpi_communicator);

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
    sum (const SparseMatrix<Number> &local,
         const MPI_Comm &mpi_communicator,
         SparseMatrix<Number> &global);

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
    T max (const T &t,
           const MPI_Comm &mpi_communicator);

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
    void max (const T        &values,
              const MPI_Comm &mpi_communicator,
              U              &maxima);

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
    void max (const ArrayView<const T> &values,
              const MPI_Comm           &mpi_communicator,
              const ArrayView<T>       &maxima);

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
    T min (const T &t,
           const MPI_Comm &mpi_communicator);

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
    void min (const T        &values,
              const MPI_Comm &mpi_communicator,
              U              &minima);

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
    void min (const ArrayView<const T> &values,
              const MPI_Comm           &mpi_communicator,
              const ArrayView<T>       &minima);

    /**
     * A data structure to store the result of the min_max_avg() function.
     * The structure stores the minimum, maximum, and average of one
     * value contributed by each processor that participates in an
     * @ref GlossMPICommunicator "MPI communicator". The structure also stores
     * the indices (or, more precisely, the @ref GlossMPIRank "MPI rank")
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
       * One of the ranks (i.e., @ref GlossMPIRank "MPI rank" within an
       * @ref GlossMPICommunicator "MPI communicator") of the
       * processors that hold the minimal value.
       */
      unsigned int min_index;

      /**
       * One of the ranks (i.e., @ref GlossMPIRank "MPI rank" within an
       * @ref GlossMPICommunicator "MPI communicator") of the
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
    min_max_avg (const double my_value,
                 const MPI_Comm &mpi_communicator);

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
      MPI_InitFinalize (int    &argc,
                        char ** &argv,
                        const unsigned int max_num_threads = numbers::invalid_unsigned_int);

      /**
       * Destructor. Calls <tt>MPI_Finalize()</tt> in case this class owns the
       * MPI process.
       */
      ~MPI_InitFinalize();
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
    bool job_supports_mpi ();

    /**
     * Initiate a some-to-some communication, and exchange arbitrary objects
     * (the class T should be serializable using boost::serialize) between
     * processors.
     *
     * @param[in] comm MPI communicator.
     *
     * @param[in] objects_to_send A map from the rank (unsigned int) of the
     *  process meant to receive the data and the object to send (the type T
     *  must be serializable for this function to work properly).
     *
     * @return A map from the rank (unsigned int) of the process
     *  which sent the data and object received.
     *
     * @author Giovanni Alzetta, Luca Heltai, 2017
     */
    template <typename T>
    std::map<unsigned int, T>
    some_to_some(const MPI_Comm                                &comm,
                 const std::map <unsigned int, T>     &objects_to_send);

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
     *
     * @author Giovanni Alzetta, Luca Heltai, 2017
     */
    template <typename T>
    std::vector<T>
    all_gather(const MPI_Comm &comm,
               const T        &object_to_send);

    /**
     * A generalization of the classic MPI_Gather function, that accepts
     * arbitrary data types T, as long as boost::serialize accepts T as an
     * argument.
     *
     * @param[in] comm MPI communicator.
     * @param[in] object_to_send an object to send to the root process
     * @param[in] root_process The process, which receives the objects from all processes.
     *  By default the process with rank 0 is the root process.
     *
     * @return The @p root_process receives a vector of objects, with size equal to the number of
     *  processes in the MPI communicator. Each entry contains the object
     *  received from the processor with the corresponding rank within the
     *  communicator. All other processes receive an empty vector.
     *
     * @author Benjamin Brands, 2017
     */
    template <typename T>
    std::vector<T>
    gather(const MPI_Comm &comm,
           const T &object_to_send,
           const unsigned int root_process=0);

#ifndef DOXYGEN
    // declaration for an internal function that lives in mpi.templates.h
    namespace internal
    {
      template <typename T>
      void all_reduce (const MPI_Op             &mpi_op,
                       const ArrayView<const T> &values,
                       const MPI_Comm           &mpi_communicator,
                       const ArrayView<T>       &output);
    }

    // Since these depend on N they must live in the header file
    template <typename T, unsigned int N>
    void sum (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&sums)[N])
    {
      internal::all_reduce(MPI_SUM, ArrayView<const T>(values, N),
                           mpi_communicator, ArrayView<T>(sums, N));
    }

    template <typename T, unsigned int N>
    void max (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&maxima)[N])
    {
      internal::all_reduce(MPI_MAX, ArrayView<const T>(values, N),
                           mpi_communicator, ArrayView<T>(maxima, N));
    }

    template <typename T, unsigned int N>
    void min (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&minima)[N])
    {
      internal::all_reduce(MPI_MIN, ArrayView<const T>(values, N),
                           mpi_communicator, ArrayView<T>(minima, N));
    }

    template<typename T>
    std::map<unsigned int, T>
    some_to_some(const MPI_Comm &comm,
                 const std::map<unsigned int, T> &objects_to_send)
    {
#ifndef DEAL_II_WITH_MPI
      (void)comm;
      Assert(objects_to_send.size() == 0, ExcMessage("Cannot send to more than one processor."));
      Assert(objects_to_send.find(0) != objects_to_send.end() || objects_to_send.size() == 0,
             ExcMessage("Can only send to myself or to nobody."));
      return objects_to_send;
#else

      std::vector<unsigned int> send_to(objects_to_send.size());
      {
        unsigned int i=0;
        for (const auto &m: objects_to_send)
          send_to[i++] = m.first;
      }
      AssertDimension(send_to.size(), objects_to_send.size());

      const auto receive_from =
        Utilities::MPI::compute_point_to_point_communication_pattern(comm, send_to);

      // Sending buffers
      std::vector<std::vector<char> > buffers_to_send(send_to.size());
      std::vector<MPI_Request> buffer_send_requests(send_to.size());
      {
        unsigned int i = 0;
        for (const auto &rank_obj : objects_to_send)
          {
            const auto &rank = rank_obj.first;
            buffers_to_send[i] = Utilities::pack(rank_obj.second);
            const int ierr = MPI_Isend(buffers_to_send[i].data(),
                                       buffers_to_send[i].size(), MPI_CHAR,
                                       rank, 21, comm, &buffer_send_requests[i]);
            AssertThrowMPI(ierr);
            ++i;
          }
      }

      // Receiving buffers
      std::map<unsigned int, T> received_objects;
      {
        std::vector<char> buffer;
        // We do this on a first come/first served basis
        for (unsigned int i = 0; i<receive_from.size(); ++i)
          {
            // Probe what's going on. Take data from the first available sender
            MPI_Status status;
            int ierr = MPI_Probe(MPI_ANY_SOURCE, 21, comm, &status);
            AssertThrowMPI(ierr);

            // Length of the message
            int len;
            ierr = MPI_Get_count(&status, MPI_CHAR, &len);
            AssertThrowMPI(ierr);
            buffer.resize(len);

            // Source rank
            const unsigned int rank = status.MPI_SOURCE;

            // Actually receive the message
            ierr = MPI_Recv(buffer.data(), len, MPI_CHAR,
                            rank, 21, comm, MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);
            Assert(received_objects.find(rank) == received_objects.end(),
                   ExcInternalError("I should not receive again from this rank"));
            received_objects[rank] = Utilities::unpack<T>(buffer);
          }
      }

      // Wait to have sent all objects.
      MPI_Waitall(send_to.size(), buffer_send_requests.data(),MPI_STATUSES_IGNORE);

      return received_objects;
#endif // deal.II with MPI
    }

    template<typename T>
    std::vector<T> all_gather(const MPI_Comm       &comm,
                              const T &object)
    {
#ifndef DEAL_II_WITH_MPI
      (void)comm;
      std::vector<T> v(1, object);
      return v;
#else
      const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);

      std::vector<char> buffer = Utilities::pack(object);

      int n_local_data = buffer.size();

      // Vector to store the size of loc_data_array for every process
      std::vector<int> size_all_data(n_procs,0);

      // Exchanging the size of each buffer
      MPI_Allgather(&n_local_data, 1, MPI_INT,
                    &(size_all_data[0]), 1, MPI_INT,
                    comm);

      // Now computing the displacement, relative to recvbuf,
      // at which to store the incoming buffer
      std::vector<int> rdispls(n_procs);
      rdispls[0] = 0;
      for (unsigned int i=1; i < n_procs; ++i)
        rdispls[i] = rdispls[i-1] + size_all_data[i-1];

      // Step 3: exchange the buffer:
      std::vector<char> received_unrolled_buffer(rdispls.back() + size_all_data.back());

      MPI_Allgatherv(buffer.data(), n_local_data, MPI_CHAR,
                     received_unrolled_buffer.data(), size_all_data.data(),
                     rdispls.data(), MPI_CHAR, comm);

      std::vector<T> received_objects(n_procs);
      for (unsigned int i= 0; i < n_procs; ++i)
        {
          std::vector<char> local_buffer(received_unrolled_buffer.begin()+rdispls[i],
                                         received_unrolled_buffer.begin()+rdispls[i]+size_all_data[i]);
          received_objects[i] = Utilities::unpack<T>(local_buffer);
        }

      return received_objects;
#endif
    }

    template <typename T>
    std::vector<T>
    gather(const MPI_Comm &comm,
           const T &object_to_send,
           const unsigned int root_process)
    {
#ifndef DEAL_II_WITH_MPI
      (void)comm;
      (void)root_process;
      std::vector<T> v(1, object_to_send);
      return v;
#else
      const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(comm);
      const auto my_rank = dealii::Utilities::MPI::this_mpi_process(comm);

      Assert(root_process < n_procs, ExcIndexRange(root_process,0,n_procs));

      std::vector<char> buffer = Utilities::pack(object_to_send);
      int n_local_data = buffer.size();

      // Vector to store the size of loc_data_array for every process
      // only the root process needs to allocate memory for that purpose
      std::vector<int> size_all_data;
      if (my_rank==root_process)
        size_all_data.resize(n_procs,0);

      // Exchanging the size of each buffer
      int ierr = MPI_Gather(&n_local_data, 1, MPI_INT,
                            size_all_data.data(), 1, MPI_INT,
                            root_process, comm);
      AssertThrowMPI(ierr);

      // Now computing the displacement, relative to recvbuf,
      // at which to store the incoming buffer; only for root
      std::vector<int> rdispls;
      if (my_rank==root_process)
        {
          rdispls.resize(n_procs,0);
          for (unsigned int i=1; i<n_procs; ++i)
            rdispls[i] = rdispls[i-1] + size_all_data[i-1];
        }
      // exchange the buffer:
      std::vector<char> received_unrolled_buffer;
      if (my_rank==root_process)
        received_unrolled_buffer.resize(rdispls.back() + size_all_data.back());

      ierr = MPI_Gatherv(buffer.data(), n_local_data, MPI_CHAR,
                         received_unrolled_buffer.data(), size_all_data.data(),
                         rdispls.data(), MPI_CHAR,
                         root_process, comm);
      AssertThrowMPI(ierr);

      std::vector<T> received_objects;

      if (my_rank==root_process)
        {
          received_objects.resize(n_procs);

          for (unsigned int i=0; i<n_procs; ++i)
            {
              const std::vector<char> local_buffer(received_unrolled_buffer.begin()+rdispls[i],
                                                   received_unrolled_buffer.begin()+rdispls[i]+size_all_data[i]);
              received_objects[i] = Utilities::unpack<T>(local_buffer);
            }
        }
      return received_objects;
#endif
    }

#endif
  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
