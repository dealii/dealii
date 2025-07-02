// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_la_parallel_vector_h
#define dealii_la_parallel_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/communication_pattern_base.h>
#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/memory_space.h>
#include <deal.II/base/memory_space_data.h>
#include <deal.II/base/mpi_stub.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/partitioner.h>

#include <deal.II/lac/read_vector.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_type_traits.h>

#include <iomanip>
#include <memory>

#if defined(DEAL_II_WITH_MPI)
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <mpi.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#endif


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
namespace LinearAlgebra
{
  /**
   * A namespace for parallel implementations of vectors.
   */
  namespace distributed
  {
    template <typename, typename>
    class BlockVector;
  }

  template <typename>
  class ReadWriteVector;
} // namespace LinearAlgebra

#  ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  namespace MPI
  {
    class Vector;
  }
} // namespace PETScWrappers
#  endif

#  ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class Vector;
  }
} // namespace TrilinosWrappers
#  endif
#endif

namespace LinearAlgebra
{
  namespace distributed
  {
    /**
     * @addtogroup Vectors
     * @{
     */

    /**
     * Implementation of a parallel vector class. The design of this class is
     * similar to the standard ::dealii::Vector class in deal.II, with the
     * exception that storage is distributed with MPI.
     *
     * The vector is designed for the following scheme of parallel
     * partitioning:
     * <ul>
     * <li> The indices held by individual processes (locally owned part) in
     * the MPI parallelization form a contiguous range
     * <code>[my_first_index,my_last_index)</code>.
     * <li> Ghost indices residing on arbitrary positions of other processors
     * are allowed. It is in general more efficient if ghost indices are
     * clustered, since they are stored as a set of intervals. The
     * communication pattern of the ghost indices is determined when calling
     * the function <code>reinit (locally_owned, ghost_indices,
     * communicator)</code>, and retained until the partitioning is changed.
     * This allows for efficient parallel communication of indices. In
     * particular, it stores the communication pattern, rather than having to
     * compute it again for every communication. For more information on ghost
     * vectors, see also the
     * @ref GlossGhostedVector "glossary entry on vectors with ghost elements".
     * <li> Besides the usual global access operator() it is also possible to
     * access vector entries in the local index space with the function @p
     * local_element(). Locally owned indices are placed first, [0,
     * locally_owned_size()), and then all ghost indices follow after them
     * contiguously, [locally_owned_size(),
     * locally_owned_size()+get_partitioner()->n_ghost_indices()).
     * </ul>
     *
     * Functions related to parallel functionality:
     * <ul>
     * <li> The function <code>compress()</code> goes through the data
     * associated with ghost indices and communicates it to the owner process,
     * which can then add it to the correct position. This can be used e.g.
     * after having run an assembly routine involving ghosts that fill this
     * vector. Note that the @p insert mode of @p compress() does not set the
     * elements included in ghost entries but simply discards them, assuming
     * that the owning processor has set them to the desired value already
     * (See also the
     * @ref GlossCompress "glossary entry on compress").
     * <li> The <code>update_ghost_values()</code> function imports the data
     * from the owning processor to the ghost indices in order to provide read
     * access to the data associated with ghosts.
     * <li> It is possible to split the above functions into two phases, where
     * the first initiates the communication and the second one finishes it.
     * These functions can be used to overlap communication with computations
     * in other parts of the code.
     * <li> Of course, reduction operations (like norms) make use of
     * collective all-to-all MPI communications.
     * </ul>
     *
     * This vector can take two different states with respect to ghost
     * elements:
     * <ul>
     * <li> After creation and whenever zero_out_ghost_values() is called (or
     * <code>operator= (0.)</code>), the vector does only allow writing into
     * ghost elements but not reading from ghost elements.
     * <li> After a call to update_ghost_values(), the vector does not allow
     * writing into ghost elements but only reading from them. This is to
     * avoid undesired ghost data artifacts when calling compress() after
     * modifying some vector entries. The current status of the ghost entries
     * (read mode or write mode) can be queried by the method
     * has_ghost_elements(), which returns <code>true</code> exactly when
     * ghost elements have been updated and <code>false</code> otherwise,
     * irrespective of the actual number of ghost entries in the vector layout
     * (for that information, use
     * <code>get_partitioner()->n_ghost_indices()</code> instead).
     * </ul>
     *
     * This vector uses the facilities of the class dealii::Vector<Number> for
     * implementing the operations on the local range of the vector. In
     * particular, it also inherits thread parallelism that splits most
     * vector-vector operations into smaller chunks if the program uses
     * multiple threads. This may or may not be desired when working also with
     * MPI.
     *
     * <h4>Limitations regarding the vector size</h4>
     *
     * This vector class is based on two different number types for indexing.
     * The so-called global index type encodes the overall size of the vector.
     * Its type is types::global_dof_index. The largest possible value is
     * <code>2^32-1</code> or approximately 4 billion in case 64 bit integers
     * are disabled at configuration of deal.II (default case) or
     * <code>2^64-1</code> or approximately <code>10^19</code> if 64 bit
     * integers are enabled (see the glossary entry on
     * @ref GlobalDoFIndex
     * for further information).
     *
     * The second relevant index type is the local index used within one MPI
     * rank. As opposed to the global index, the implementation assumes 32-bit
     * unsigned integers unconditionally. In other words, to actually use a
     * vector with more than four billion entries, you need to use MPI with
     * more than one rank (which in general is a safe assumption since four
     * billion entries consume at least 16 GB of memory for floats or 32 GB of
     * memory for doubles) and enable 64-bit indices. If more than 4 billion
     * local elements are present, the implementation tries to detect that,
     * which triggers an exception and aborts the code. Note, however, that
     * the detection of overflow is tricky and the detection mechanism might
     * fail in some circumstances. Therefore, it is strongly recommended to
     * not rely on this class to automatically detect the unsupported case.
     *
     * <h4>GPU support</h4>
     *
     * This vector class supports two different memory spaces: Host and Default.
     * If the MemorySpace template argument is not specified, the memory space
     * is Host and all the data is allocated on the CPU. When the memory space
     * is Default, all the data is allocated on Kokkos' default memory space.
     * That means that if Kokkos was configured with a GPU backend, the data is
     * allocated on a GPU. The operations on the vector are performed on the
     * chosen memory space. From the host, there are two methods to access the
     * elements of the Vector when using the Default memory space:
     * <ul>
     * <li> use get_values():
     * @code
     * Vector<double, MemorySpace::Default> vector(local_range, comm);
     * double* vector_dev = vector.get_values();
     * const int n_local_elements = local_range.n_elements();
     * std::vector<double> vector_host(n_local_elements, 1.);
     * Kokkos::deep_copy(Kokkos::View<double, Kokkos::HostSpace>(
     *                     vector_host.data(), n_local_elements),
     *                   Kokkos::View<double,
     * MemorySpace::Default::kokkos_space>( vector_dev, n_local_elements));
     * @endcode
     * <li> use import_elements():
     * @code
     * Vector<double, MemorySpace::Default> vector(local_range, comm);
     * ReadWriteVector<double> rw_vector(local_range);
     * for (auto & val : rw_vector)
     *   val = 1.;
     * vector.import_elements(rw_vector, VectorOperations::insert);
     * @endcode
     * </ul>
     * The import method is a lot safer and will perform an MPI communication if
     * necessary. Since an MPI communication may be performed, import needs to
     * be called on all the processors.
     *
     * @note By default, the GPU @ref GlossDevice "device" id is chosen in a round-robin fashion
     * according to the local MPI rank id. To choose a different @ref GlossDevice "device", Kokkos
     * has to be initialized explicitly providing the respective device id
     * explicitly.
     *
     * <h4>MPI-3 shared-memory support</h4>
     *
     * In Host mode, this class allows to use MPI-3 shared-memory features
     * by providing a separate MPI communicator that consists of processes on
     * the same shared-memory domain. By calling
     * `vector.shared_vector_data();`,
     * users have read-only access to both locally-owned and ghost values of
     * processes combined in the shared-memory communicator (@p comm_sm in
     * reinit()).
     *
     * For this to work, you have to call the constructor or one of the reinit()
     * functions of this class with a non-default value for the `comm_sm`
     * argument, where the argument corresponds to a communicator consisting of
     * all processes on the same shared-memory domain. This kind of communicator
     * can be created using the following code snippet:
     * @code
     *   MPI_Comm comm_sm;
     *   MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
     *                       &comm_sm);
     * @endcode
     */
    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class Vector : public ::dealii::ReadVector<Number>
    {
    public:
      using memory_space    = MemorySpace;
      using value_type      = Number;
      using pointer         = value_type *;
      using const_pointer   = const value_type *;
      using iterator        = value_type *;
      using const_iterator  = const value_type *;
      using reference       = value_type &;
      using const_reference = const value_type &;
      using size_type       = types::global_dof_index;
      using real_type       = typename numbers::NumberTraits<Number>::real_type;

      static_assert(
        std::is_same_v<MemorySpace, ::dealii::MemorySpace::Host> ||
          std::is_same_v<MemorySpace, ::dealii::MemorySpace::Default>,
        "MemorySpace should be Host or Default");

      static_assert(
        (!std::is_same_v<MemorySpace, ::dealii::MemorySpace::Default>) ||
          std::is_same_v<Number, float> || std::is_same_v<Number, double>,
        "Number should be float or double for Default memory space");

      /**
       * @name 1: Basic Object-handling
       */
      /** @{ */
      /**
       * Empty constructor.
       */
      Vector();

      /**
       * Copy constructor. Uses the parallel partitioning of @p in_vector.
       * It should be noted that this constructor automatically sets ghost
       * values to zero. Call @p update_ghost_values() directly following
       * construction if a ghosted vector is required.
       */
      Vector(const Vector<Number, MemorySpace> &in_vector);

      /**
       * Move constructor. Uses the swap method.
       *
       * @note In order for this constructor to leave the moved-from object in a
       * valid state it must allocate memory (in this case, an empty
       * partitioner) - hence it cannot be marked as noexcept.
       */
      Vector(Vector<Number, MemorySpace> &&in_vector); // NOLINT

      /**
       * Construct a parallel vector of the given global size without any
       * actual parallel distribution.
       */
      Vector(const size_type size);

      /**
       * Construct a parallel vector. The local range is specified by @p
       * locally_owned_set (note that this must be a contiguous interval,
       * multiple intervals are not possible). The IndexSet @p ghost_indices
       * specifies ghost indices, i.e., indices which one might need to read
       * data from or accumulate data from. It is allowed that the set of
       * ghost indices also contains the local range, but it does not need to.
       *
       * This function involves global communication, so it should only be
       * called once for a given layout. Use the constructor with
       * Vector<Number> argument to create additional vectors with the same
       * parallel layout.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      Vector(const IndexSet &local_range,
             const IndexSet &ghost_indices,
             const MPI_Comm  communicator);

      /**
       * Same constructor as above but without any ghost indices.
       */
      Vector(const IndexSet &local_range, const MPI_Comm communicator);

      /**
       * Create the vector based on the parallel partitioning described in @p
       * partitioner. The input argument is a shared pointer, which stores the
       * partitioner data only once and share it between several vectors with
       * the same layout.
       */
      Vector(
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

      /**
       * Destructor.
       */
      ~Vector();

      /**
       * Set the global size of the vector to @p size without any actual
       * parallel distribution.
       */
      void
      reinit(const size_type size, const bool omit_zeroing_entries = false);

      /**
       * Uses the parallel layout of the input vector @p in_vector and
       * allocates memory for this vector. Recommended initialization function
       * when several vectors with the same layout should be created.
       *
       * If the flag @p omit_zeroing_entries is set to false, the memory will
       * be initialized with zero, otherwise the memory will be untouched (and
       * the user must make sure to fill it with reasonable data before using
       * it).
       */
      template <typename Number2>
      void
      reinit(const Vector<Number2, MemorySpace> &in_vector,
             const bool                          omit_zeroing_entries = false);

      /**
       * Initialize the vector. The local range is specified by @p local_range
       * (note that this must be a contiguous interval, multiple intervals are
       * not possible). The IndexSet @p ghost_indices specifies ghost indices,
       * i.e., indices which one might need to read data from or accumulate data
       * from. It is allowed that the set of ghost indices also contains the
       * local range, but it does not need to.
       *
       * This function involves global communication, so it should only be
       * called once for a given layout. Use the @p reinit function with
       * Vector<Number> argument to create additional vectors with the same
       * parallel layout.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      void
      reinit(const IndexSet &local_range,
             const IndexSet &ghost_indices,
             const MPI_Comm  communicator);

      /**
       * Same as above, but without ghost entries.
       */
      void
      reinit(const IndexSet &local_range, const MPI_Comm communicator);

      /**
       * Initialize the vector given to the parallel partitioning described in
       * @p partitioner. The input argument is a shared pointer, which stores
       * the partitioner data only once and can be shared between several
       * vectors with the same layout.
       *
       * The optional argument @p comm_sm, which consists of processes on
       * the same shared-memory domain, allows users have read-only access to
       * both locally-owned and ghost values of processes combined in the
       * shared-memory communicator. See the general documentation of this class
       * for more information about this argument.
       */
      void
      reinit(
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
        const MPI_Comm comm_sm = MPI_COMM_SELF);

      /**
       * This function exists purely for reasons of compatibility with the
       * PETScWrappers::MPI::Vector and TrilinosWrappers::MPI::Vector classes.
       *
       * It calls the function above, and ignores the parameter @p make_ghosted.
       */
      void
      reinit(
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
        const bool                                                make_ghosted,
        const MPI_Comm &comm_sm = MPI_COMM_SELF);

      /**
       * Initialize vector with @p local_size locally-owned and @p ghost_size
       * ghost degrees of freedoms.
       *
       * The optional argument @p comm_sm, which consists of processes on
       * the same shared-memory domain, allows users have read-only access to
       * both locally-owned and ghost values of processes combined in the
       * shared-memory communicator. See the general documentation of this class
       * for more information about this argument.
       *
       * @note In the created underlying partitioner, the local index range is
       *   translated to global indices in an ascending and one-to-one fashion,
       *   i.e., the indices of process $p$ sit exactly between the indices of
       *   the processes $p-1$ and $p+1$, respectively. Setting the
       *   @p ghost_size variable to an appropriate value provides memory space
       *   for the ghost data in a vector's memory allocation as and allows
       *   access to it via local_element(). However, the associated global
       *   indices must be handled externally in this case.
       */
      void
      reinit(const types::global_dof_index local_size,
             const types::global_dof_index ghost_size,
             const MPI_Comm                comm,
             const MPI_Comm                comm_sm = MPI_COMM_SELF);

      /**
       * Swap the contents of this vector and the other vector @p v. One could
       * do this operation with a temporary variable and copying over the data
       * elements, but this function is significantly more efficient since it
       * only swaps the pointers to the data of the two vectors and therefore
       * does not need to allocate temporary storage and move data around.
       *
       * This function is analogous to the @p swap function of all C++
       * standard containers. Also, there is a global function
       * <tt>swap(u,v)</tt> that simply calls <tt>u.swap(v)</tt>, again in
       * analogy to standard functions.
       */
      void
      swap(Vector<Number, MemorySpace> &v) noexcept;

      /**
       * Move assignment operator.
       *
       * @note This method may throw an exception (should an MPI check fail) and
       * is consequently not `noexcept`.
       */
      Vector<Number, MemorySpace> &
      operator=(Vector<Number, MemorySpace> &&in_vector); // NOLINT

      /**
       * Assigns the vector to the parallel partitioning of the input vector
       * @p in_vector, and copies all the data.
       *
       * The semantics of this operator are complex. If the two vectors have
       * the same size, and
       * if either the left or right hand side vector of the assignment (i.e.,
       * either the input vector on the right hand side, or the calling vector
       * to the left of the assignment operator) currently has ghost elements,
       * then the left hand side vector will also have ghost values and will
       * consequently be a read-only vector (see also the
       * @ref GlossGhostedVector "glossary entry" on the issue). Otherwise, the
       * left hand vector will be a writable vector after this operation.
       * These semantics facilitate having a vector with ghost elements on the
       * left hand side of the assignment, and a vector without ghost elements
       * on the right hand side, with the resulting left hand side vector
       * having the correct values in both its locally owned and its ghost
       * elements.
       *
       * On the other hand, if the left hand side vector does not have the
       * correct size yet, or is perhaps an entirely uninitialized vector,
       * then the assignment is simply a copy operation in the usual sense:
       * In that case, if the right hand side has no ghost elements (i.e.,
       * is a completely distributed vector), then the left hand side will
       * have no ghost elements either. And if the right hand side has
       * ghost elements (and is consequently read-only), then the left
       * hand side will have these same properties after the operation.
       */
      Vector<Number, MemorySpace> &
      operator=(const Vector<Number, MemorySpace> &in_vector);

      /**
       * Assigns the vector to the parallel partitioning of the input vector
       * @p in_vector, and copies all the data.
       *
       * If one of the input vector or the calling vector (to the left of the
       * assignment operator) had ghost elements set before this operation,
       * the calling vector will have ghost values set. Otherwise, it will be
       * in write mode. If the input vector does not have any ghost elements
       * at all, the vector will also update its ghost values in analogy to
       * the respective setting the Trilinos and PETSc vectors.
       */
      template <typename Number2>
      Vector<Number, MemorySpace> &
      operator=(const Vector<Number2, MemorySpace> &in_vector);

      /** @} */

      /**
       * @name 2: Parallel data exchange
       */
      /** @{ */
      /**
       * This function copies the data that has accumulated in the data buffer
       * for ghost indices to the owning processor. For the meaning of the
       * argument @p operation, see the entry on
       * @ref GlossCompress "Compressing distributed vectors and matrices"
       * in the glossary.
       *
       * There are four variants for this function. If called with argument @p
       * VectorOperation::add adds all the data accumulated in ghost elements
       * to the respective elements on the owning processor and clears the
       * ghost array afterwards. If called with argument @p
       * VectorOperation::insert, a set operation is performed. Since setting
       * elements in a vector with ghost elements is ambiguous (as one can set
       * both the element on the ghost site as well as the owning site), this
       * operation makes the assumption that all data is set correctly on the
       * owning processor. Upon call of compress(VectorOperation::insert), all
       * ghost entries are thus simply zeroed out (using zero_ghost_values()).
       * In debug mode, a check is performed for whether the data set is
       * actually consistent between processors, i.e., whenever a non-zero
       * ghost element is found, it is compared to the value on the owning
       * processor and an exception is thrown if these elements do not agree.
       * If called with VectorOperation::min or VectorOperation::max, the
       * minimum or maximum on all elements across the processors is set.
       * @note This vector class has a fixed set of ghost entries attached to
       * the local representation. As a consequence, all ghost entries are
       * assumed to be valid and will be exchanged unconditionally according
       * to the given VectorOperation. Make sure to initialize all ghost
       * entries with the neutral element of the given VectorOperation or
       * touch all ghost entries. The neutral element is zero for
       * VectorOperation::add and VectorOperation::insert, `+inf` for
       * VectorOperation::min, and `-inf` for VectorOperation::max. If all
       * values are initialized with values below zero and compress is called
       * with VectorOperation::max two times subsequently, the maximal value
       * after the second calculation will be zero.
       */
      void
      compress(VectorOperation::values operation);

      /**
       * Fills the data field for ghost indices with the values stored in the
       * respective positions of the owning processor. This function is needed
       * before reading from ghosts. The function is @p const even though
       * ghost data is changed. This is needed to allow functions with a @p
       * const vector to perform the data exchange without creating
       * temporaries.
       *
       * After calling this method, write access to ghost elements of the
       * vector is forbidden and an exception is thrown. Only read access to
       * ghost elements is allowed in this state. Note that all subsequent
       * operations on this vector, like global vector addition, etc., will
       * also update the ghost values by a call to this method after the
       * operation. However, global reduction operations like norms or the
       * inner product will always ignore ghost elements in order to avoid
       * counting the ghost data more than once. To allow writing to ghost
       * elements again, call zero_out_ghost_values().
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      void
      update_ghost_values() const;

      /**
       * Initiates communication for the @p compress() function with
       * non-blocking communication. This function does not wait for the
       * transfer to finish, in order to allow for other computations during the
       * time it takes until all data arrives.
       *
       * Before the data is actually exchanged, the function must be followed
       * by a call to @p compress_finish().
       *
       * In case this function is called for more than one vector before @p
       * compress_finish() is invoked, it is mandatory to specify a unique
       * communication channel to each such call, in order to avoid several
       * messages with the same ID that will corrupt this operation. Any
       * communication channel less than 100 is a valid value (in particular,
       * the range $[100, 200)$ is reserved for
       * LinearAlgebra::distributed::BlockVector).
       */
      void
      compress_start(const unsigned int      communication_channel = 0,
                     VectorOperation::values operation = VectorOperation::add);

      /**
       * For all requests that have been initiated in compress_start, wait for
       * the communication to finish. Once it is finished, add or set the data
       * (depending on the flag operation) to the respective positions in the
       * owning processor, and clear the contents in the ghost data fields.
       * The meaning of this argument is the same as in compress().
       *
       * This function should be called exactly once per vector after calling
       * compress_start, otherwise the result is undefined. In particular, it
       * is not well-defined to call compress_start on the same vector again
       * before compress_finished has been called. However, there is no
       * warning to prevent this situation.
       *
       * Must follow a call to the @p compress_start function.
       *
       * When the MemorySpace is Default and MPI is not GPU-aware, data changed
       * on the @ref GlossDevice "device" after the call to compress_start will be lost.
       */
      void
      compress_finish(VectorOperation::values operation);

      /**
       * Initiates communication for the @p update_ghost_values() function
       * with non-blocking communication. This function does not wait for the
       * transfer to finish, in order to allow for other computations during
       * the time it takes until all data arrives.
       *
       * Before the data is actually exchanged, the function must be followed
       * by a call to @p update_ghost_values_finish().
       *
       * In case this function is called for more than one vector before @p
       * update_ghost_values_finish() is invoked, it is mandatory to specify a
       * unique communication channel to each such call, in order to avoid
       * several messages with the same ID that will corrupt this operation.
       * Any communication channel less than 100 is a valid value (in
       * particular, the range $[100, 200)$ is reserved for
       * LinearAlgebra::distributed::BlockVector).
       */
      void
      update_ghost_values_start(
        const unsigned int communication_channel = 0) const;


      /**
       * For all requests that have been started in update_ghost_values_start,
       * wait for the communication to finish.
       *
       * Must follow a call to the @p update_ghost_values_start function
       * before reading data from ghost indices.
       */
      void
      update_ghost_values_finish() const;

      /**
       * This method zeros the entries on ghost dofs, but does not touch
       * locally owned DoFs.
       *
       * After calling this method, read access to ghost elements of the
       * vector is forbidden and an exception is thrown. Only write access to
       * ghost elements is allowed in this state.
       */
      void
      zero_out_ghost_values() const;

      /**
       * Return whether the vector currently is in a state where ghost values
       * can be read or not. This is the same functionality as other parallel
       * vectors have. If this method returns false, this only means that
       * read-access to ghost elements is prohibited whereas write access is
       * still possible (to those entries specified as ghosts during
       * initialization), not that there are no ghost elements at all.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      bool
      has_ghost_elements() const;

      /**
       * This method copies the data in the locally owned range from another
       * distributed vector @p src into the calling vector. As opposed to
       * operator= that also includes ghost entries, this operation ignores
       * the ghost range. The only prerequisite is that the local range on the
       * calling vector and the given vector @p src are the same on all
       * processors. It is explicitly allowed that the two vectors have
       * different ghost elements that might or might not be related to each
       * other.
       *
       * Since no data exchange is performed, make sure that neither @p src
       * nor the calling vector have pending communications in order to obtain
       * correct results.
       */
      template <typename Number2>
      void
      copy_locally_owned_data_from(const Vector<Number2, MemorySpace> &src);

      /**
       * Import all the elements present in the distributed vector @p src.
       * VectorOperation::values @p operation is used to decide if the elements
       * in @p V should be added to the current vector or replace the current
       * elements. The main purpose of this function is to get data from one
       * memory space, e.g. Default, to the other, e.g. the Host.
       *
       * @note The partitioners of the two distributed vectors need to be the
       * same as no MPI communication is performed.
       */
      template <typename MemorySpace2>
      void
      import_elements(const Vector<Number, MemorySpace2> &src,
                      VectorOperation::values             operation);

      /**
       * @deprecated Use import_elements() instead.
       */
      template <typename MemorySpace2>
      DEAL_II_DEPRECATED void
      import(const Vector<Number, MemorySpace2> &src,
             VectorOperation::values             operation)
      {
        import_elements(src, operation);
      }

      /** @} */

      /**
       * @name 3: Implementation of vector space operations
       */
      /** @{ */

      /**
       * Multiply the entire vector by a fixed factor.
       */
      Vector<Number, MemorySpace> &
      operator*=(const Number factor);

      /**
       * Divide the entire vector by a fixed factor.
       */
      Vector<Number, MemorySpace> &
      operator/=(const Number factor);

      /**
       * Add the vector @p V to the present one.
       */
      Vector<Number, MemorySpace> &
      operator+=(const Vector<Number, MemorySpace> &V);

      /**
       * Subtract the vector @p V from the present one.
       */
      Vector<Number, MemorySpace> &
      operator-=(const Vector<Number, MemorySpace> &V);

      /**
       * Import all the elements present in the vector's IndexSet from the input
       * vector @p V. VectorOperation::values @p operation is used to decide if
       * the elements in @p V should be added to the current vector or replace the
       * current elements. The last parameter can be used if the same
       * communication pattern is used multiple times. This can be used to
       * improve performance.
       *
       * @note If the MemorySpace is Default, the data in the ReadWriteVector will
       * be moved to the @ref GlossDevice "device".
       */
      void
      import_elements(
        const LinearAlgebra::ReadWriteVector<Number> &V,
        const VectorOperation::values                 operation,
        const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
          &communication_pattern = {});

      /**
       * @deprecated Use import_elements() instead.
       */
      DEAL_II_DEPRECATED void
      import(const LinearAlgebra::ReadWriteVector<Number> &V,
             VectorOperation::values                       operation,
             std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
               communication_pattern = {})
      {
        import_elements(V, operation, communication_pattern);
      }

      /**
       * Return the scalar product of two vectors.
       */
      Number
      operator*(const Vector<Number, MemorySpace> &V) const;

      /**
       * Add @p a to all components. Note that @p a is a scalar not a vector.
       */
      void
      add(const Number a);

      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
       */
      void
      add(const Number a, const Vector<Number, MemorySpace> &V);

      /**
       * Multiple addition of scaled vectors, i.e. <tt>*this += a*V+b*W</tt>.
       */
      void
      add(const Number                       a,
          const Vector<Number, MemorySpace> &V,
          const Number                       b,
          const Vector<Number, MemorySpace> &W);

      /**
       * A collective add operation: This function adds a whole set of values
       * stored in @p values to the vector components specified by @p indices.
       */
      void
      add(const std::vector<size_type> &indices,
          const std::vector<Number>    &values);

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this =
       * s*(*this)+a*V</tt>.
       */
      void
      sadd(const Number                       s,
           const Number                       a,
           const Vector<Number, MemorySpace> &V);

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument. This function is mostly meant to simulate multiplication (and
       * immediate re-assignment) by a diagonal scaling matrix.
       */
      void
      scale(const Vector<Number, MemorySpace> &scaling_factors);

      /**
       * Assignment <tt>*this = a*V</tt>.
       */
      void
      equ(const Number a, const Vector<Number, MemorySpace> &V);

      /**
       * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
       * absolute values of all entries among all processors).
       */
      real_type
      l1_norm() const;

      /**
       * Return the $l_2$ norm of the vector (i.e., the square root of
       * the sum of the square of all entries among all processors).
       */
      real_type
      l2_norm() const;

      /**
       * Return the square of the $l_2$ norm of the vector.
       */
      real_type
      norm_sqr() const;

      /**
       * Return the maximum norm of the vector (i.e., the maximum absolute value
       * among all entries and among all processors).
       */
      real_type
      linfty_norm() const;

      /**
       * Perform a combined operation of a vector addition and a subsequent
       * inner product, returning the value of the inner product. In other
       * words, the result of this function is the same as if the user called
       * @code
       * this->add(a, V);
       * return_value = *this * W;
       * @endcode
       *
       * The reason this function exists is that this operation involves less
       * memory transfer than calling the two functions separately. This method
       * only needs to load three vectors, @p this, @p V, @p W, whereas calling
       * separate methods means to load the calling vector @p this twice. Since
       * most vector operations are memory transfer limited, this reduces the
       * time by 25\% (or 50\% if @p W equals @p this).
       *
       * For complex-valued vectors, the scalar product in the second step is
       * implemented as
       * $\left<v,w\right>=\sum_i v_i \bar{w_i}$.
       */
      Number
      add_and_dot(const Number                       a,
                  const Vector<Number, MemorySpace> &V,
                  const Vector<Number, MemorySpace> &W);

      /**
       * Return the global size of the vector, equal to the sum of the number of
       * locally owned indices among all processors.
       */
      virtual size_type
      size() const override;

      /**
       * Return an index set that describes which elements of this vector are
       * owned by the current processor. As a consequence, the index sets
       * returned on different processors if this is a distributed vector will
       * form disjoint sets that add up to the complete index set. Obviously, if
       * a vector is created on only one processor, then the result would
       * satisfy
       * @code
       *  vec.locally_owned_elements() == complete_index_set(vec.size())
       * @endcode
       */
      dealii::IndexSet
      locally_owned_elements() const;

      /**
       * Print the vector to the output stream @p out.
       */
      void
      print(std::ostream      &out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const;

      /**
       * Return the memory consumption of this class in bytes.
       */
      std::size_t
      memory_consumption() const;
      /** @} */

      /**
       * @name 4: Other vector operations
       */
      /** @{ */

      /**
       * Sets all elements of the vector to the scalar @p s. If the scalar is
       * zero, also ghost elements are set to zero, otherwise they remain
       * unchanged.
       */
      Vector<Number, MemorySpace> &
      operator=(const Number s);

      /**
       * This is a collective add operation that adds a whole set of values
       * stored in @p values to the vector components specified by @p indices.
       */
      template <typename OtherNumber>
      void
      add(const std::vector<size_type>        &indices,
          const ::dealii::Vector<OtherNumber> &values);

      /**
       * Take an address where n_elements are stored contiguously and add them
       * into the vector.
       */
      template <typename OtherNumber>
      void
      add(const size_type    n_elements,
          const size_type   *indices,
          const OtherNumber *values);

      /**
       * Scaling and simple vector addition, i.e.  <tt>*this =
       * s*(*this)+V</tt>.
       */
      void
      sadd(const Number s, const Vector<Number, MemorySpace> &V);

      /** @} */


      /**
       * @name 5: Entry access and local data representation
       */
      /** @{ */

      /**
       * Return the local size of the vector, i.e., the number of indices
       * owned locally.
       */
      size_type
      locally_owned_size() const;

      /**
       * Return true if the given global index is in the local range of this
       * processor.
       */
      bool
      in_local_range(const size_type global_index) const;

      /**
       * Make the @p Vector class a bit like the <tt>vector<></tt> class of
       * the C++ standard library by returning iterators to the start and end
       * of the <i>locally owned</i> elements of this vector.
       *
       * It holds that end() - begin() == locally_owned_size().
       *
       * @note For the Default memory space, the iterator might point to memory
       * on the @ref GlossDevice "device".
       */
      iterator
      begin();

      /**
       * Return constant iterator to the start of the locally owned elements
       * of the vector.
       *
       * @note For the Default memory space, the iterator might point to memory
       * on the @ref GlossDevice "device".
       */
      const_iterator
      begin() const;

      /**
       * Return an iterator pointing to the element past the end of the array
       * of locally owned entries.
       *
       * @note For the Default memory space, the iterator might point to memory
       * on the @ref GlossDevice "device".
       */
      iterator
      end();

      /**
       * Return a constant iterator pointing to the element past the end of
       * the array of the locally owned entries.
       *
       * @note For the Default memory space, the iterator might point to memory
       * on the @ref GlossDevice "device".
       */
      const_iterator
      end() const;

      /**
       * Read access to the data in the position corresponding to @p
       * global_index. The index must be either in the local range of the
       * vector or be specified as a ghost index at construction.
       *
       * Performance: <tt>O(1)</tt> for locally owned elements that represent
       * a contiguous range and <tt>O(log(n<sub>ranges</sub>))</tt> for ghost
       * elements (quite fast, but slower than local_element()).
       */
      Number
      operator()(const size_type global_index) const;

      /**
       * Read and write access to the data in the position corresponding to @p
       * global_index. The index must be either in the local range of the
       * vector or be specified as a ghost index at construction.
       *
       * Performance: <tt>O(1)</tt> for locally owned elements that represent
       * a contiguous range and <tt>O(log(n<sub>ranges</sub>))</tt> for ghost
       * elements (quite fast, but slower than local_element()).
       */
      Number &
      operator()(const size_type global_index);

      /**
       * Read access to the data in the position corresponding to @p
       * global_index. The index must be either in the local range of the
       * vector or be specified as a ghost index at construction.
       *
       * This function does the same thing as operator().
       */
      Number
      operator[](const size_type global_index) const;
      /**
       * Read and write access to the data in the position corresponding to @p
       * global_index. The index must be either in the local range of the
       * vector or be specified as a ghost index at construction.
       *
       * This function does the same thing as operator().
       */
      Number &
      operator[](const size_type global_index);

      /**
       * Read access to the data field specified by @p local_index. Locally
       * owned indices can be accessed with indices
       * <code>[0,locally_owned_size)</code>, and ghost indices with indices
       * <code>[locally_owned_size,locally_owned_size+get_partitioner()->n_ghost_indices()]</code>.
       *
       * Performance: Direct array access (fast).
       */
      Number
      local_element(const size_type local_index) const;

      /**
       * Read and write access to the data field specified by @p local_index.
       * Locally owned indices can be accessed with indices
       * <code>[0,locally_owned_size())</code>, and ghost indices with indices
       * <code>[locally_owned_size(), locally_owned_size()+n_ghosts]</code>.
       *
       * Performance: Direct array access (fast).
       */
      Number &
      local_element(const size_type local_index);

      /**
       * Return the pointer to the underlying raw array.
       *
       * @note For the Default memory space, the pointer might point to memory
       * on the @ref GlossDevice "device".
       */
      Number *
      get_values() const;

      /**
       * Instead of getting individual elements of a vector via operator(),
       * this function allows getting a whole set of elements at once. The
       * indices of the elements to be read are stated in the first argument,
       * the corresponding values are returned in the second.
       *
       * If the current vector is called @p v, then this function is the equivalent
       * to the code
       * @code
       *   for (unsigned int i=0; i<indices.size(); ++i)
       *     values[i] = v[indices[i]];
       * @endcode
       *
       * @pre The sizes of the @p indices and @p values arrays must be identical.
       *
       * @note This function is not implemented for Default memory space.
       */
      template <typename OtherNumber>
      void
      extract_subvector_to(const std::vector<size_type> &indices,
                           std::vector<OtherNumber>     &values) const;

      /**
       * Extract a range of elements all at once.
       */
      virtual void
      extract_subvector_to(
        const ArrayView<const types::global_dof_index> &indices,
        const ArrayView<Number> &elements) const override;

      /**
       * Instead of getting individual elements of a vector via operator(),
       * this function allows getting a whole set of elements at once. In
       * contrast to the previous function, this function obtains the
       * indices of the elements by dereferencing all elements of the iterator
       * range provided by the first two arguments, and puts the vector
       * values into memory locations obtained by dereferencing a range
       * of iterators starting at the location pointed to by the third
       * argument.
       *
       * If the current vector is called @p v, then this function is the equivalent
       * to the code
       * @code
       *   ForwardIterator indices_p = indices_begin;
       *   OutputIterator  values_p  = values_begin;
       *   while (indices_p != indices_end)
       *   {
       *     *values_p = v[*indices_p];
       *     ++indices_p;
       *     ++values_p;
       *   }
       * @endcode
       *
       * @pre It must be possible to write into as many memory locations
       *   starting at @p values_begin as there are iterators between
       *   @p indices_begin and @p indices_end.
       */
      template <typename ForwardIterator, typename OutputIterator>
      void
      extract_subvector_to(ForwardIterator       indices_begin,
                           const ForwardIterator indices_end,
                           OutputIterator        values_begin) const;
      /**
       * Return whether the vector contains only elements with value zero.
       * This is a @ref GlossCollectiveOperation "collective operation". This function is expensive, because
       * potentially all elements have to be checked.
       */
      bool
      all_zero() const;

      /**
       * Compute the mean value of all the entries in the vector.
       */
      Number
      mean_value() const;

      /**
       * $l_p$-norm of the vector. The pth root of the sum of the pth powers
       * of the absolute values of the elements.
       */
      real_type
      lp_norm(const real_type p) const;
      /** @} */

      /**
       * @name 6: Mixed stuff
       */
      /** @{ */

      /**
       * Return the underlying MPI communicator.
       */
      MPI_Comm
      get_mpi_communicator() const;

      /**
       * Return the MPI partitioner that describes the parallel layout of the
       * vector. This object can be used to initialize another vector with the
       * respective reinit() call, for additional queries regarding the
       * parallel communication, or the compatibility of partitioners.
       */
      const std::shared_ptr<const Utilities::MPI::Partitioner> &
      get_partitioner() const;

      /**
       * Check whether the given partitioner is compatible with the
       * partitioner used for this vector. Two partitioners are compatible if
       * they have the same local size and the same ghost indices. They do not
       * necessarily need to be the same data field of the shared pointer.
       * This is a local operation only, i.e., if only some processors decide
       * that the partitioning is not compatible, only these processors will
       * return @p false, whereas the other processors will return @p true.
       */
      bool
      partitioners_are_compatible(
        const Utilities::MPI::Partitioner &part) const;

      /**
       * Check whether the given partitioner is compatible with the
       * partitioner used for this vector. Two partitioners are compatible if
       * they have the same local size and the same ghost indices. They do not
       * necessarily need to be the same data field. As opposed to
       * partitioners_are_compatible(), this method checks for compatibility
       * among all processors and the method only returns @p true if the
       * partitioner is the same on all processors.
       *
       * This method performs global communication, so make sure to use it
       * only in a context where all processors call it the same number of
       * times.
       */
      bool
      partitioners_are_globally_compatible(
        const Utilities::MPI::Partitioner &part) const;

      /**
       * Change the ghost state of this vector to @p ghosted.
       */
      void
      set_ghost_state(const bool ghosted) const;

      /**
       * Get pointers to the beginning of the values of the other
       * processes of the same shared-memory domain.
       */
      const std::vector<ArrayView<const Number>> &
      shared_vector_data() const;

      /** @} */

      /**
       * Attempt to perform an operation between two incompatible vector types.
       *
       * @ingroup Exceptions
       */
      DeclException0(ExcVectorTypeNotCompatible);

      /**
       * Exception
       */
      DeclException3(ExcNonMatchingElements,
                     Number,
                     Number,
                     unsigned int,
                     << "Called compress(VectorOperation::insert), but"
                     << " the element received from a remote processor, value "
                     << std::setprecision(16) << arg1
                     << ", does not match with the value "
                     << std::setprecision(16) << arg2
                     << " on the owner processor " << arg3);

      /**
       * Exception
       */
      DeclException4(
        ExcAccessToNonLocalElement,
        size_type,
        size_type,
        size_type,
        size_type,
        << "You tried to access element " << arg1
        << " of a distributed vector, but this element is not "
        << "stored on the current processor. Note: The range of "
        << "locally owned elements is [" << arg2 << ',' << arg3
        << "], and there are " << arg4 << " ghost elements "
        << "that this vector can access."
        << "\n\n"
        << "A common source for this kind of problem is that you "
        << "are passing a 'fully distributed' vector into a function "
        << "that needs read access to vector elements that correspond "
        << "to degrees of freedom on ghost cells (or at least to "
        << "'locally active' degrees of freedom that are not also "
        << "'locally owned'). You need to pass a vector that has these "
        << "elements as ghost entries.");

    private:
      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>
       * without MPI communication.
       */
      void
      add_local(const Number a, const Vector<Number, MemorySpace> &V);

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this =
       * s*(*this)+a*V</tt> without MPI communication.
       */
      void
      sadd_local(const Number                       s,
                 const Number                       a,
                 const Vector<Number, MemorySpace> &V);

      /**
       * Local part of the inner product of two vectors.
       */
      template <typename Number2>
      Number
      inner_product_local(const Vector<Number2, MemorySpace> &V) const;

      /**
       * Local part of norm_sqr().
       */
      real_type
      norm_sqr_local() const;

      /**
       * Local part of mean_value().
       */
      Number
      mean_value_local() const;

      /**
       * Local part of l1_norm().
       */
      real_type
      l1_norm_local() const;

      /**
       * Local part of lp_norm().
       */
      real_type
      lp_norm_local(const real_type p) const;

      /**
       * Local part of linfty_norm().
       */
      real_type
      linfty_norm_local() const;

      /**
       * Local part of the addition followed by an inner product of two
       * vectors. The same applies for complex-valued vectors as for
       * the add_and_dot() function.
       */
      Number
      add_and_dot_local(const Number                       a,
                        const Vector<Number, MemorySpace> &V,
                        const Vector<Number, MemorySpace> &W);

      /**
       * Assert that there are no spurious non-zero entries in the ghost
       * region of the vector caused by some forgotten compress() or
       * zero_out_ghost_values() calls, which is an invariant of the vector
       * space operations such as the addition of vectors, scaling a vector,
       * and similar.
       */
      void
      assert_no_residual_content_in_ghost_region() const;

      /**
       * Shared pointer to store the parallel partitioning information. This
       * information can be shared between several vectors that have the same
       * partitioning.
       */
      std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;

      /**
       * The size that is currently allocated in the val array.
       */
      size_type allocated_size;

      /**
       * Underlying data structure storing the local elements of this vector.
       */
      mutable ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> data;

      /**
       * For parallel loops with TBB, this member variable stores the affinity
       * information of loops.
       */
      mutable std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
        thread_loop_partitioner;

      /**
       * Temporary storage that holds the data that is sent to this processor
       * in compress() or sent from this processor in update_ghost_values().
       */
      mutable ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace>
        import_data;

      /**
       * Stores whether the vector currently allows for reading ghost elements
       * or not. Note that this is to ensure consistent ghost data and does
       * not indicate whether the vector actually can store ghost elements. In
       * particular, when assembling a vector we do not allow reading
       * elements, only writing them.
       */
      mutable bool vector_is_ghosted;

#ifdef DEAL_II_WITH_MPI
      /**
       * A vector that collects all requests from compress() operations.
       * This class uses persistent MPI communicators, i.e., the communication
       * channels are stored during successive calls to a given function. This
       * reduces the overhead involved with setting up the MPI machinery, but
       * it does not remove the need for a receive operation to be posted
       * before the data can actually be sent.
       */
      std::vector<MPI_Request> compress_requests;

      /**
       * A vector that collects all requests from update_ghost_values()
       * operations. This class uses persistent MPI communicators.
       */
      mutable std::vector<MPI_Request> update_ghost_values_requests;
#endif

      /**
       * A lock that makes sure that the compress() and update_ghost_values()
       * functions give reasonable results also when used
       * with several threads.
       */
      mutable std::mutex mutex;

      /**
       * Communicator to be used for the shared-memory domain. See the general
       * documentation of this class for more information about the purpose of
       * `comm_sm`.
       */
      MPI_Comm comm_sm;

      /**
       * A helper function that clears the compress_requests and
       * update_ghost_values_requests field. Used in reinit() functions.
       */
      void
      clear_mpi_requests();

      /**
       * A helper function that is used to resize the val array.
       */
      void
      resize_val(const size_type new_allocated_size,
                 const MPI_Comm  comm_sm = MPI_COMM_SELF);

      // Make all other vector types friends.
      template <typename Number2, typename MemorySpace2>
      friend class Vector;

      // Make BlockVector type friends.
      template <typename Number2, typename MemorySpace2>
      friend class BlockVector;
    };
    /** @} */


    /*-------------------- Inline functions ---------------------------------*/

#ifndef DOXYGEN

    template <typename Number, typename MemorySpace>
    inline bool
    Vector<Number, MemorySpace>::has_ghost_elements() const
    {
      return vector_is_ghosted;
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::size_type
    Vector<Number, MemorySpace>::size() const
    {
      return partitioner->size();
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::size_type
    Vector<Number, MemorySpace>::locally_owned_size() const
    {
      return partitioner->locally_owned_size();
    }



    template <typename Number, typename MemorySpace>
    inline bool
    Vector<Number, MemorySpace>::in_local_range(
      const size_type global_index) const
    {
      return partitioner->in_local_range(global_index);
    }



    template <typename Number, typename MemorySpace>
    inline IndexSet
    Vector<Number, MemorySpace>::locally_owned_elements() const
    {
      IndexSet is(size());

      is.add_range(partitioner->local_range().first,
                   partitioner->local_range().second);

      return is;
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::iterator
    Vector<Number, MemorySpace>::begin()
    {
      return data.values.data();
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::const_iterator
    Vector<Number, MemorySpace>::begin() const
    {
      return data.values.data();
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::iterator
    Vector<Number, MemorySpace>::end()
    {
      return data.values.data() + partitioner->locally_owned_size();
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::const_iterator
    Vector<Number, MemorySpace>::end() const
    {
      return data.values.data() + partitioner->locally_owned_size();
    }



    template <typename Number, typename MemorySpace>
    const std::vector<ArrayView<const Number>> &
    Vector<Number, MemorySpace>::shared_vector_data() const
    {
      return data.values_sm;
    }



    template <typename Number, typename MemorySpace>
    inline Number
    Vector<Number, MemorySpace>::operator()(const size_type global_index) const
    {
      Assert((std::is_same_v<MemorySpace, ::dealii::MemorySpace::Host>),
             ExcMessage(
               "This function is only implemented for the Host memory space"));
      Assert(
        partitioner->in_local_range(global_index) ||
          partitioner->ghost_indices().is_element(global_index),
        ExcAccessToNonLocalElement(global_index,
                                   partitioner->local_range().first,
                                   partitioner->local_range().second == 0 ?
                                     0 :
                                     (partitioner->local_range().second - 1),
                                   partitioner->ghost_indices().n_elements()));
      // do not allow reading a vector which is not in ghost mode
      Assert(partitioner->in_local_range(global_index) ||
               vector_is_ghosted == true,
             ExcMessage("You tried to read a ghost element of this vector, "
                        "but it has not imported its ghost values."));
      return data.values[partitioner->global_to_local(global_index)];
    }



    template <typename Number, typename MemorySpace>
    inline Number &
    Vector<Number, MemorySpace>::operator()(const size_type global_index)
    {
      Assert((std::is_same_v<MemorySpace, ::dealii::MemorySpace::Host>),
             ExcMessage(
               "This function is only implemented for the Host memory space"));
      Assert(
        partitioner->in_local_range(global_index) ||
          partitioner->ghost_indices().is_element(global_index),
        ExcAccessToNonLocalElement(global_index,
                                   partitioner->local_range().first,
                                   partitioner->local_range().second == 0 ?
                                     0 :
                                     (partitioner->local_range().second - 1),
                                   partitioner->ghost_indices().n_elements()));
      // we would like to prevent reading ghosts from a vector that does not
      // have them imported, but this is not possible because we might be in a
      // part of the code where the vector has enabled ghosts but is non-const
      // (then, the compiler picks this method according to the C++ rule book
      // even if a human would pick the const method when this subsequent use
      // is just a read)
      return data.values[partitioner->global_to_local(global_index)];
    }



    template <typename Number, typename MemorySpace>
    inline Number
    Vector<Number, MemorySpace>::operator[](const size_type global_index) const
    {
      return operator()(global_index);
    }



    template <typename Number, typename MemorySpace>
    inline Number &
    Vector<Number, MemorySpace>::operator[](const size_type global_index)
    {
      return operator()(global_index);
    }



    template <typename Number, typename MemorySpace>
    inline Number
    Vector<Number, MemorySpace>::local_element(
      const size_type local_index) const
    {
      Assert((std::is_same_v<MemorySpace, ::dealii::MemorySpace::Host>),
             ExcMessage(
               "This function is only implemented for the Host memory space"));
      AssertIndexRange(local_index,
                       partitioner->locally_owned_size() +
                         partitioner->n_ghost_indices());
      // do not allow reading a vector which is not in ghost mode
      Assert(local_index < locally_owned_size() || vector_is_ghosted == true,
             ExcMessage("You tried to read a ghost element of this vector, "
                        "but it has not imported its ghost values."));

      return data.values[local_index];
    }



    template <typename Number, typename MemorySpace>
    inline Number &
    Vector<Number, MemorySpace>::local_element(const size_type local_index)
    {
      Assert((std::is_same_v<MemorySpace, ::dealii::MemorySpace::Host>),
             ExcMessage(
               "This function is only implemented for the Host memory space"));

      AssertIndexRange(local_index,
                       partitioner->locally_owned_size() +
                         partitioner->n_ghost_indices());

      return data.values[local_index];
    }



    template <typename Number, typename MemorySpace>
    inline Number *
    Vector<Number, MemorySpace>::get_values() const
    {
      return data.values.data();
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::extract_subvector_to(
      const std::vector<size_type> &indices,
      std::vector<OtherNumber>     &values) const
    {
      for (size_type i = 0; i < indices.size(); ++i)
        values[i] = operator()(indices[i]);
    }



    template <typename Number, typename MemorySpace>
    template <typename ForwardIterator, typename OutputIterator>
    inline void
    Vector<Number, MemorySpace>::extract_subvector_to(
      ForwardIterator       indices_begin,
      const ForwardIterator indices_end,
      OutputIterator        values_begin) const
    {
      while (indices_begin != indices_end)
        {
          *values_begin = operator()(*indices_begin);
          ++indices_begin;
          ++values_begin;
        }
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::add(
      const std::vector<size_type>        &indices,
      const ::dealii::Vector<OtherNumber> &values)
    {
      AssertDimension(indices.size(), values.size());
      for (size_type i = 0; i < indices.size(); ++i)
        {
          Assert(
            numbers::is_finite(values[i]),
            ExcMessage(
              "The given value is not finite but either infinite or Not A Number (NaN)"));
          this->operator()(indices[i]) += values(i);
        }
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::add(const size_type    n_elements,
                                     const size_type   *indices,
                                     const OtherNumber *values)
    {
      for (size_type i = 0; i < n_elements; ++i, ++indices, ++values)
        {
          Assert(
            numbers::is_finite(*values),
            ExcMessage(
              "The given value is not finite but either infinite or Not A Number (NaN)"));
          this->operator()(*indices) += *values;
        }
    }



    template <typename Number, typename MemorySpace>
    inline MPI_Comm
    Vector<Number, MemorySpace>::get_mpi_communicator() const
    {
      return partitioner->get_mpi_communicator();
    }



    template <typename Number, typename MemorySpace>
    inline const std::shared_ptr<const Utilities::MPI::Partitioner> &
    Vector<Number, MemorySpace>::get_partitioner() const
    {
      return partitioner;
    }



    template <typename Number, typename MemorySpace>
    inline void
    Vector<Number, MemorySpace>::set_ghost_state(const bool ghosted) const
    {
      vector_is_ghosted = ghosted;
    }

#endif

  } // namespace distributed
} // namespace LinearAlgebra


/**
 * Global function @p swap which overloads the default implementation of the
 * C++ standard library which uses a temporary object. The function simply
 * exchanges the data of the two vectors.
 *
 * @relatesalso Vector
 */
template <typename Number, typename MemorySpace>
inline void
swap(LinearAlgebra::distributed::Vector<Number, MemorySpace> &u,
     LinearAlgebra::distributed::Vector<Number, MemorySpace> &v) noexcept
{
  u.swap(v);
}


/**
 * Declare dealii::LinearAlgebra::distributed::Vector as distributed vector.
 */
template <typename Number, typename MemorySpace>
struct is_serial_vector<LinearAlgebra::distributed::Vector<Number, MemorySpace>>
  : std::false_type
{};



namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename>
    class ReinitHelper;

    /**
     * A helper class used internally in linear_operator.h. Specialization for
     * LinearAlgebra::distributed::Vector.
     */
    template <typename Number>
    class ReinitHelper<LinearAlgebra::distributed::Vector<Number>>
    {
    public:
      // A helper type-trait that leverage SFINAE to figure out if type T has
      // void T::get_mpi_communicator()
      template <typename T>
      using get_mpi_communicator_t =
        decltype(std::declval<T>().get_mpi_communicator());

      template <typename T>
      static constexpr bool has_get_mpi_communicator =
        is_supported_operation<get_mpi_communicator_t, T>;

      // A helper type-trait that leverage SFINAE to figure out if type T has
      // void T::locally_owned_domain_indices()
      template <typename T>
      using locally_owned_domain_indices_t =
        decltype(std::declval<T>().locally_owned_domain_indices());

      template <typename T>
      static constexpr bool has_locally_owned_domain_indices =
        is_supported_operation<locally_owned_domain_indices_t, T>;

      // A helper type-trait that leverage SFINAE to figure out if type T has
      // void T::locally_owned_range_indices()
      template <typename T>
      using locally_owned_range_indices_t =
        decltype(std::declval<T>().locally_owned_range_indices());

      template <typename T>
      static constexpr bool has_locally_owned_range_indices =
        is_supported_operation<locally_owned_range_indices_t, T>;

      // A helper type-trait that leverage SFINAE to figure out if type T has
      // void T::initialize_dof_vector(VectorType v)
      template <typename T>
      using initialize_dof_vector_t =
        decltype(std::declval<T>().initialize_dof_vector(
          std::declval<LinearAlgebra::distributed::Vector<Number> &>()));

      template <typename T>
      static constexpr bool has_initialize_dof_vector =
        is_supported_operation<initialize_dof_vector_t, T>;

      // Used for (Trilinos/PETSc)Wrappers::SparseMatrix
      template <typename MatrixType,
#if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1900
                std::enable_if_t<has_get_mpi_communicator<MatrixType> &&
                                   has_locally_owned_domain_indices<MatrixType>,
#else
                // workaround for Intel 18
                std::enable_if_t<
                  is_supported_operation<get_mpi_communicator_t, MatrixType> &&
                    is_supported_operation<locally_owned_domain_indices_t,
                                           MatrixType>,
#endif
                                 MatrixType> * = nullptr>
      static void
      reinit_domain_vector(MatrixType                                 &mat,
                           LinearAlgebra::distributed::Vector<Number> &vec,
                           bool /*omit_zeroing_entries*/)
      {
        vec.reinit(mat.locally_owned_domain_indices(),
                   mat.get_mpi_communicator());
      }

      // Used for MatrixFree and DiagonalMatrix
      template <typename MatrixType,
#if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1900
                std::enable_if_t<has_initialize_dof_vector<MatrixType>,
#else
                // workaround for Intel 18
                std::enable_if_t<
                  is_supported_operation<initialize_dof_vector_t, MatrixType>,
#endif
                                 MatrixType> * = nullptr>
      static void
      reinit_domain_vector(MatrixType                                 &mat,
                           LinearAlgebra::distributed::Vector<Number> &vec,
                           bool omit_zeroing_entries)
      {
        mat.initialize_dof_vector(vec);
        if (!omit_zeroing_entries)
          vec = Number();
      }

      // Used for (Trilinos/PETSc)Wrappers::SparseMatrix
      template <typename MatrixType,
#if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1900
                std::enable_if_t<has_get_mpi_communicator<MatrixType> &&
                                   has_locally_owned_range_indices<MatrixType>,
#else
                // workaround for Intel 18
                std::enable_if_t<
                  is_supported_operation<get_mpi_communicator_t, MatrixType> &&
                    is_supported_operation<locally_owned_range_indices_t,
                                           MatrixType>,
#endif
                                 MatrixType> * = nullptr>
      static void
      reinit_range_vector(MatrixType                                 &mat,
                          LinearAlgebra::distributed::Vector<Number> &vec,
                          bool /*omit_zeroing_entries*/)
      {
        vec.reinit(mat.locally_owned_range_indices(),
                   mat.get_mpi_communicator());
      }

      // Used for MatrixFree and DiagonalMatrix
      template <typename MatrixType,
#if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1900
                std::enable_if_t<has_initialize_dof_vector<MatrixType>,
#else
                // workaround for Intel 18
                std::enable_if_t<
                  is_supported_operation<initialize_dof_vector_t, MatrixType>,
#endif
                                 MatrixType> * = nullptr>
      static void
      reinit_range_vector(MatrixType                                 &mat,
                          LinearAlgebra::distributed::Vector<Number> &vec,
                          bool omit_zeroing_entries)
      {
        mat.initialize_dof_vector(vec);
        if (!omit_zeroing_entries)
          vec = Number();
      }
    };

  } // namespace LinearOperatorImplementation
} /* namespace internal */


DEAL_II_NAMESPACE_CLOSE

#endif
