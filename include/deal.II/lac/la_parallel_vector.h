// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2019 by the deal.II authors
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

#ifndef dealii_la_parallel_vector_h
#define dealii_la_parallel_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_space_vector.h>
#include <deal.II/lac/vector_type_traits.h>

#include <iomanip>
#include <memory>

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
    template <typename>
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
    /*! @addtogroup Vectors
     *@{
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
     * local_size()), and then all ghost indices follow after them
     * contiguously, [local_size(), local_size()+n_ghost_entries()).
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
     * <li> After creation and whenever zero_out_ghosts() is called (or
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
     * (for that information, use n_ghost_entries() instead).
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
     * <h4>CUDA support</h4>
     *
     * This vector class supports two different memory spaces: Host and CUDA. By
     * default, the memory space is Host and all the data are allocated on the
     * CPU. When the memory space is CUDA, all the data is allocated on the GPU.
     * The operations on the vector are performed on the chosen memory space. *
     * From the host, there are two methods to access the elements of the Vector
     * when using the CUDA memory space:
     * <ul>
     * <li> use get_values():
     * @code
     * Vector<double, MemorySpace::CUDA> vector(local_range, comm);
     * double* vector_dev = vector.get_values();
     * std::vector<double> vector_host(local_range.n_elements(), 1.);
     * Utilities::CUDA::copy_to_dev(vector_host, vector_dev);
     * @endcode
     * <li> use import():
     * @code
     * Vector<double, MemorySpace::CUDA> vector(local_range, comm);
     * ReadWriteVector<double> rw_vector(local_range);
     * for (auto & val : rw_vector)
     *   val = 1.;
     * vector.import(rw_vector, VectorOperations::insert);
     * @endcode
     * </ul>
     * The import method is a lot safer and will perform an MPI communication if
     * necessary. Since an MPI communication may be performed, import needs to
     * be called on all the processors.
     *
     * @note By default, all the ranks will try to access the device 0. This is
     * fine is if you have one rank per node and one gpu per node. If you
     * have multiple GPUs on one node, we need each process to access a
     * different GPU. If each node has the same number of GPUs, this can be done
     * as follows:
     * <code> int n_devices = 0; cudaGetDeviceCount(&n_devices); int
     * device_id = my_rank % n_devices;
     * cudaSetDevice(device_id);
     * </code>
     * @see CUDAWrappers
     *
     * @author Katharina Kormann, Martin Kronbichler, Bruno Turcksin 2010, 2011,
     * 2016, 2018
     */
    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class Vector : public ::dealii::LinearAlgebra::VectorSpaceVector<Number>,
                   public Subscriptor
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
        std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value ||
          std::is_same<MemorySpace, ::dealii::MemorySpace::CUDA>::value,
        "MemorySpace should be Host or CUDA");

      /**
       * @name 1: Basic Object-handling
       */
      //@{
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
       * partitioner. The input argument is a shared pointer, which store the
       * partitioner data only once and share it between several vectors with
       * the same layout.
       */
      Vector(
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

      /**
       * Destructor.
       */
      virtual ~Vector() override;

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
       * Initialize the vector. The local range is specified by @p
       * locally_owned_set (note that this must be a contiguous interval,
       * multiple intervals are not possible). The IndexSet @p ghost_indices
       * specifies ghost indices, i.e., indices which one might need to read
       * data from or accumulate data from. It is allowed that the set of
       * ghost indices also contains the local range, but it does not need to.
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
       * @p partitioner. The input argument is a shared pointer, which store
       * the partitioner data only once and share it between several vectors
       * with the same layout.
       */
      void
      reinit(
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

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
       *
       * This function is virtual in order to allow for derived classes to
       * handle memory separately.
       */
      void
      swap(Vector<Number, MemorySpace> &v);

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

      //@}

      /**
       * @name 2: Parallel data exchange
       */
      //@{
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
      virtual void
      compress(::dealii::VectorOperation::values operation) override;

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
       * elements again, call zero_out_ghosts().
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      void
      update_ghost_values() const;

      /**
       * Initiates communication for the @p compress() function with non-
       * blocking communication. This function does not wait for the transfer
       * to finish, in order to allow for other computations during the time
       * it takes until all data arrives.
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
      compress_start(
        const unsigned int                communication_channel = 0,
        ::dealii::VectorOperation::values operation = VectorOperation::add);

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
       */
      void
      compress_finish(::dealii::VectorOperation::values operation);

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
      zero_out_ghosts() const;

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
       * memory space, e.g. CUDA, to the other, e.g. the Host.
       *
       * @note The partitioners of the two distributed vectors need to be the
       * same as no MPI communication is performed.
       */
      template <typename MemorySpace2>
      void
      import(const Vector<Number, MemorySpace2> &src,
             VectorOperation::values             operation);

      //@}

      /**
       * @name 3: Implementation of VectorSpaceVector
       */
      //@{

      /**
       * Change the dimension to that of the vector V. The elements of V are not
       * copied.
       */
      virtual void
      reinit(const VectorSpaceVector<Number> &V,
             const bool omit_zeroing_entries = false) override;

      /**
       * Multiply the entire vector by a fixed factor.
       */
      virtual Vector<Number, MemorySpace> &
      operator*=(const Number factor) override;

      /**
       * Divide the entire vector by a fixed factor.
       */
      virtual Vector<Number, MemorySpace> &
      operator/=(const Number factor) override;

      /**
       * Add the vector @p V to the present one.
       */
      virtual Vector<Number, MemorySpace> &
      operator+=(const VectorSpaceVector<Number> &V) override;

      /**
       * Subtract the vector @p V from the present one.
       */
      virtual Vector<Number, MemorySpace> &
      operator-=(const VectorSpaceVector<Number> &V) override;

      /**
       * Import all the elements present in the vector's IndexSet from the input
       * vector @p V. VectorOperation::values @p operation is used to decide if
       * the elements in @p V should be added to the current vector or replace the
       * current elements. The last parameter can be used if the same
       * communication pattern is used multiple times. This can be used to
       * improve performance.
       *
       * @note If the MemorySpace is CUDA, the data in the ReadWriteVector will
       * be moved to the device.
       */
      virtual void
      import(
        const LinearAlgebra::ReadWriteVector<Number> &  V,
        VectorOperation::values                         operation,
        std::shared_ptr<const CommunicationPatternBase> communication_pattern =
          std::shared_ptr<const CommunicationPatternBase>()) override;

      /**
       * Return the scalar product of two vectors.
       */
      virtual Number
      operator*(const VectorSpaceVector<Number> &V) const override;

      /**
       * Add @p a to all components. Note that @p a is a scalar not a vector.
       */
      virtual void
      add(const Number a) override;

      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
       */
      virtual void
      add(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * Multiple addition of scaled vectors, i.e. <tt>*this += a*V+b*W</tt>.
       */
      virtual void
      add(const Number                     a,
          const VectorSpaceVector<Number> &V,
          const Number                     b,
          const VectorSpaceVector<Number> &W) override;

      /**
       * A collective add operation: This function adds a whole set of values
       * stored in @p values to the vector components specified by @p indices.
       */
      virtual void
      add(const std::vector<size_type> &indices,
          const std::vector<Number> &   values);

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this =
       * s*(*this)+a*V</tt>.
       */
      virtual void
      sadd(const Number                     s,
           const Number                     a,
           const VectorSpaceVector<Number> &V) override;

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument. This function is mostly meant to simulate multiplication (and
       * immediate re-assignment) by a diagonal scaling matrix.
       */
      virtual void
      scale(const VectorSpaceVector<Number> &scaling_factors) override;

      /**
       * Assignment <tt>*this = a*V</tt>.
       */
      virtual void
      equ(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
       * absolute values of all entries among all processors).
       */
      virtual real_type
      l1_norm() const override;

      /**
       * Return the $l_2$ norm of the vector (i.e., the square root of
       * the sum of the square of all entries among all processors).
       */
      virtual real_type
      l2_norm() const override;

      /**
       * Return the square of the $l_2$ norm of the vector.
       */
      real_type
      norm_sqr() const;

      /**
       * Return the maximum norm of the vector (i.e., the maximum absolute value
       * among all entries and among all processors).
       */
      virtual real_type
      linfty_norm() const override;

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
      virtual Number
      add_and_dot(const Number                     a,
                  const VectorSpaceVector<Number> &V,
                  const VectorSpaceVector<Number> &W) override;

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
      virtual dealii::IndexSet
      locally_owned_elements() const override;

      /**
       * Print the vector to the output stream @p out.
       */
      virtual void
      print(std::ostream &     out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const override;

      /**
       * Return the memory consumption of this class in bytes.
       */
      virtual std::size_t
      memory_consumption() const override;
      //@}

      /**
       * @name 4: Other vector operations not included in VectorSpaceVector
       */
      //@{

      /**
       * Sets all elements of the vector to the scalar @p s. If the scalar is
       * zero, also ghost elements are set to zero, otherwise they remain
       * unchanged.
       */
      virtual Vector<Number, MemorySpace> &
      operator=(const Number s) override;

      /**
       * This is a collective add operation that adds a whole set of values
       * stored in @p values to the vector components specified by @p indices.
       */
      template <typename OtherNumber>
      void
      add(const std::vector<size_type> &       indices,
          const ::dealii::Vector<OtherNumber> &values);

      /**
       * Take an address where n_elements are stored contiguously and add them
       * into the vector.
       */
      template <typename OtherNumber>
      void
      add(const size_type    n_elements,
          const size_type *  indices,
          const OtherNumber *values);

      /**
       * Scaling and simple vector addition, i.e.  <tt>*this =
       * s*(*this)+V</tt>.
       */
      void
      sadd(const Number s, const Vector<Number, MemorySpace> &V);

      //@}


      /**
       * @name 5: Entry access and local data representation
       */
      //@{

      /**
       * Return the local size of the vector, i.e., the number of indices
       * owned locally.
       */
      size_type
      local_size() const;

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
       * It holds that end() - begin() == local_size().
       *
       * @note For the CUDA memory space, the iterator points to memory on the
       * device.
       */
      iterator
      begin();

      /**
       * Return constant iterator to the start of the locally owned elements
       * of the vector.
       *
       * @note For the CUDA memory space, the iterator points to memory on the
       * device.
       */
      const_iterator
      begin() const;

      /**
       * Return an iterator pointing to the element past the end of the array
       * of locally owned entries.
       *
       * @note For the CUDA memory space, the iterator points to memory on the
       * device.
       */
      iterator
      end();

      /**
       * Return a constant iterator pointing to the element past the end of
       * the array of the locally owned entries.
       *
       * @note For the CUDA memory space, the iterator points to memory on the
       * device.
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
      Number operator[](const size_type global_index) const;
      /**
       * Read and write access to the data in the position corresponding to @p
       * global_index. The index must be either in the local range of the
       * vector or be specified as a ghost index at construction.
       *
       * This function does the same thing as operator().
       */
      Number &operator[](const size_type global_index);

      /**
       * Read access to the data field specified by @p local_index. Locally
       * owned indices can be accessed with indices
       * <code>[0,local_size)</code>, and ghost indices with indices
       * <code>[local_size,local_size+ n_ghost_entries]</code>.
       *
       * Performance: Direct array access (fast).
       */
      Number
      local_element(const size_type local_index) const;

      /**
       * Read and write access to the data field specified by @p local_index.
       * Locally owned indices can be accessed with indices
       * <code>[0,local_size)</code>, and ghost indices with indices
       * <code>[local_size,local_size+n_ghosts]</code>.
       *
       * Performance: Direct array access (fast).
       */
      Number &
      local_element(const size_type local_index);

      /**
       * Return the pointer to the underlying raw array.
       *
       * @note For the CUDA memory space, the pointer points to memory on the
       * device.
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
       * @note This function is not implemented for CUDA memory space.
       */
      template <typename OtherNumber>
      void
      extract_subvector_to(const std::vector<size_type> &indices,
                           std::vector<OtherNumber> &    values) const;

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
       * This is a collective operation. This function is expensive, because
       * potentially all elements have to be checked.
       */
      virtual bool
      all_zero() const override;

      /**
       * Compute the mean value of all the entries in the vector.
       */
      virtual Number
      mean_value() const override;

      /**
       * $l_p$-norm of the vector. The pth root of the sum of the pth powers
       * of the absolute values of the elements.
       */
      real_type
      lp_norm(const real_type p) const;
      //@}

      /**
       * @name 6: Mixed stuff
       */
      //@{

      /**
       * Return a reference to the MPI communicator object in use with this
       * vector.
       */
      const MPI_Comm &
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

      //@}

      /**
       * Attempt to perform an operation between two incompatible vector types.
       *
       * @ingroup Exceptions
       */
      DeclException0(ExcVectorTypeNotCompatible);

      /**
       * Attempt to perform an operation not implemented on the device.
       *
       * @ingroup Exceptions
       */
      DeclException0(ExcNotAllowedForCuda);

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
      DeclException4(ExcAccessToNonLocalElement,
                     size_type,
                     size_type,
                     size_type,
                     size_type,
                     << "You tried to access element " << arg1
                     << " of a distributed vector, but this element is not "
                     << "stored on the current processor. Note: The range of "
                     << "locally owned elements is " << arg2 << " to " << arg3
                     << ", and there are " << arg4 << " ghost elements "
                     << "that this vector can access.");

    private:
      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>
       * without MPI communication.
       */
      void
      add_local(const Number a, const VectorSpaceVector<Number> &V);

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this =
       * s*(*this)+a*V</tt> without MPI communication.
       */
      void
      sadd_local(const Number                     s,
                 const Number                     a,
                 const VectorSpaceVector<Number> &V);

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
       * in @p compress() or sent from this processor in
       * @p update_ghost_values.
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
       * A vector that collects all requests from @p compress() operations.
       * This class uses persistent MPI communicators, i.e., the communication
       * channels are stored during successive calls to a given function. This
       * reduces the overhead involved with setting up the MPI machinery, but
       * it does not remove the need for a receive operation to be posted
       * before the data can actually be sent.
       */
      std::vector<MPI_Request> compress_requests;

      /**
       * A vector that collects all requests from @p update_ghost_values()
       * operations. This class uses persistent MPI communicators.
       */
      mutable std::vector<MPI_Request> update_ghost_values_requests;
#endif

      /**
       * A lock that makes sure that the @p compress and @p
       * update_ghost_values functions give reasonable results also when used
       * with several threads.
       */
      mutable std::mutex mutex;

      /**
       * A helper function that clears the compress_requests and
       * update_ghost_values_requests field. Used in reinit functions.
       */
      void
      clear_mpi_requests();

      /**
       * A helper function that is used to resize the val array.
       */
      void
      resize_val(const size_type new_allocated_size);

      // Make all other vector types friends.
      template <typename Number2, typename MemorySpace2>
      friend class Vector;

      // Make BlockVector type friends.
      template <typename Number2>
      friend class BlockVector;
    };
    /*@}*/


    /*-------------------- Inline functions ---------------------------------*/

#ifndef DOXYGEN

    namespace internal
    {
      template <typename Number, typename MemorySpace>
      struct Policy
      {
        static inline typename Vector<Number, MemorySpace>::iterator
        begin(::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> &)
        {
          return nullptr;
        }

        static inline typename Vector<Number, MemorySpace>::const_iterator
        begin(
          const ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> &)
        {
          return nullptr;
        }

        static inline Number *
        get_values(
          ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace> &)
        {
          return nullptr;
        }
      };



      template <typename Number>
      struct Policy<Number, ::dealii::MemorySpace::Host>
      {
        static inline
          typename Vector<Number, ::dealii::MemorySpace::Host>::iterator
          begin(::dealii::MemorySpace::
                  MemorySpaceData<Number, ::dealii::MemorySpace::Host> &data)
        {
          return data.values.get();
        }

        static inline
          typename Vector<Number, ::dealii::MemorySpace::Host>::const_iterator
          begin(const ::dealii::MemorySpace::
                  MemorySpaceData<Number, ::dealii::MemorySpace::Host> &data)
        {
          return data.values.get();
        }

        static inline Number *
        get_values(::dealii::MemorySpace::
                     MemorySpaceData<Number, ::dealii::MemorySpace::Host> &data)
        {
          return data.values.get();
        }
      };



      template <typename Number>
      struct Policy<Number, ::dealii::MemorySpace::CUDA>
      {
        static inline
          typename Vector<Number, ::dealii::MemorySpace::CUDA>::iterator
          begin(::dealii::MemorySpace::
                  MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &data)
        {
          return data.values_dev.get();
        }

        static inline
          typename Vector<Number, ::dealii::MemorySpace::CUDA>::const_iterator
          begin(const ::dealii::MemorySpace::
                  MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &data)
        {
          return data.values_dev.get();
        }

        static inline Number *
        get_values(::dealii::MemorySpace::
                     MemorySpaceData<Number, ::dealii::MemorySpace::CUDA> &data)
        {
          return data.values_dev.get();
        }
      };
    } // namespace internal


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
    Vector<Number, MemorySpace>::local_size() const
    {
      return partitioner->local_size();
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
      return internal::Policy<Number, MemorySpace>::begin(data);
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::const_iterator
    Vector<Number, MemorySpace>::begin() const
    {
      return internal::Policy<Number, MemorySpace>::begin(data);
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::iterator
    Vector<Number, MemorySpace>::end()
    {
      return internal::Policy<Number, MemorySpace>::begin(data) +
             partitioner->local_size();
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::const_iterator
    Vector<Number, MemorySpace>::end() const
    {
      return internal::Policy<Number, MemorySpace>::begin(data) +
             partitioner->local_size();
    }



    template <typename Number, typename MemorySpace>
    inline Number
    Vector<Number, MemorySpace>::operator()(const size_type global_index) const
    {
      Assert((std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value),
             ExcMessage(
               "This function is only implemented for the Host memory space"));
      Assert(
        partitioner->in_local_range(global_index) ||
          partitioner->ghost_indices().is_element(global_index),
        ExcAccessToNonLocalElement(global_index,
                                   partitioner->local_range().first,
                                   partitioner->local_range().second,
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
      Assert((std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value),
             ExcMessage(
               "This function is only implemented for the Host memory space"));
      Assert(
        partitioner->in_local_range(global_index) ||
          partitioner->ghost_indices().is_element(global_index),
        ExcAccessToNonLocalElement(global_index,
                                   partitioner->local_range().first,
                                   partitioner->local_range().second,
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
    inline Number Vector<Number, MemorySpace>::
                  operator[](const size_type global_index) const
    {
      return operator()(global_index);
    }



    template <typename Number, typename MemorySpace>
    inline Number &Vector<Number, MemorySpace>::
                   operator[](const size_type global_index)
    {
      return operator()(global_index);
    }



    template <typename Number, typename MemorySpace>
    inline Number
    Vector<Number, MemorySpace>::local_element(
      const size_type local_index) const
    {
      Assert((std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value),
             ExcMessage(
               "This function is only implemented for the Host memory space"));
      AssertIndexRange(local_index,
                       partitioner->local_size() +
                         partitioner->n_ghost_indices());
      // do not allow reading a vector which is not in ghost mode
      Assert(local_index < local_size() || vector_is_ghosted == true,
             ExcMessage("You tried to read a ghost element of this vector, "
                        "but it has not imported its ghost values."));

      return data.values[local_index];
    }



    template <typename Number, typename MemorySpace>
    inline Number &
    Vector<Number, MemorySpace>::local_element(const size_type local_index)
    {
      Assert((std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value),
             ExcMessage(
               "This function is only implemented for the Host memory space"));

      AssertIndexRange(local_index,
                       partitioner->local_size() +
                         partitioner->n_ghost_indices());

      return data.values[local_index];
    }



    template <typename Number, typename MemorySpace>
    inline Number *
    Vector<Number, MemorySpace>::get_values() const
    {
      return internal::Policy<Number, MemorySpace>::get_values(data);
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::extract_subvector_to(
      const std::vector<size_type> &indices,
      std::vector<OtherNumber> &    values) const
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
          indices_begin++;
          values_begin++;
        }
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::add(
      const std::vector<size_type> &       indices,
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
                                     const size_type *  indices,
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
    inline const MPI_Comm &
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
 * @author Katharina Kormann, Martin Kronbichler, 2011
 */
template <typename Number, typename MemorySpace>
inline void
swap(LinearAlgebra::distributed::Vector<Number, MemorySpace> &u,
     LinearAlgebra::distributed::Vector<Number, MemorySpace> &v)
{
  u.swap(v);
}


/**
 * Declare dealii::LinearAlgebra::Vector as distributed vector.
 *
 * @author Uwe Koecher, 2017
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
      template <typename Matrix>
      static void
      reinit_range_vector(const Matrix &                              matrix,
                          LinearAlgebra::distributed::Vector<Number> &v,
                          bool omit_zeroing_entries)
      {
        matrix.initialize_dof_vector(v);
        if (!omit_zeroing_entries)
          v = Number();
      }

      template <typename Matrix>
      static void
      reinit_domain_vector(const Matrix &                              matrix,
                           LinearAlgebra::distributed::Vector<Number> &v,
                           bool omit_zeroing_entries)
      {
        matrix.initialize_dof_vector(v);
        if (!omit_zeroing_entries)
          v = Number();
      }
    };

  } // namespace LinearOperatorImplementation
} /* namespace internal */


DEAL_II_NAMESPACE_CLOSE

#endif
