// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_vector_h
#define dealii_trilinos_tpetra_vector_h


#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#include <deal.II/lac/trilinos_tpetra_types.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/mpi_stub.h>

#  include <deal.II/lac/read_vector.h>
#  include <deal.II/lac/trilinos_tpetra_communication_pattern.h>
#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/vector_operation.h>
#  include <deal.II/lac/vector_type_traits.h>

#  include <Teuchos_Comm.hpp>
#  include <Teuchos_OrdinalTraits.hpp>
#  include <Tpetra_Core.hpp>
#  include <Tpetra_Vector.hpp>
#  include <Tpetra_Version.hpp>

#  include <memory>
#  include <optional>

DEAL_II_NAMESPACE_OPEN

/**
 * Type trait indicating if a certain number type has been explicitly
 * instantiated in Tpetra. deal.II only supports those number types in Tpetra
 * wrapper classes.
 */
template <typename Number>
struct is_tpetra_type : std::false_type
{};

#  ifdef HAVE_TPETRA_INST_FLOAT
template <>
struct is_tpetra_type<float> : std::true_type
{};
#  endif

#  ifdef HAVE_TPETRA_INST_DOUBLE
template <>
struct is_tpetra_type<double> : std::true_type
{};
#  endif

#  ifdef DEAL_II_WITH_COMPLEX_VALUES
#    ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
template <>
struct is_tpetra_type<std::complex<float>> : std::true_type
{};
#    endif

#    ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
template <>
struct is_tpetra_type<std::complex<double>> : std::true_type
{};
#    endif
#  endif

namespace LinearAlgebra
{
  // Forward declaration
#  ifndef DOXYGEN
  template <typename Number>
  class ReadWriteVector;
#  endif

  /**
   * @addtogroup TpetraWrappers
   * @{
   */

  /**
   * A namespace for classes that provide wrappers for Trilinos' Tpetra vectors.
   *
   * This namespace provides wrappers for the Tpetra::Vector class from the
   * Tpetra package (https://trilinos.github.io/tpetra.html) that is part of
   * Trilinos.
   *
   * @ingroup TpetraWrappers
   */
  namespace TpetraWrappers
  {

    /**
     * This class defines type aliases that are used in vector classes
     * within the TpetraWrappers namespace.
     */
    class VectorTraits
    {
    public:
      using size_type = types::global_dof_index;
    };

    /**
     * @cond internal
     */

    /**
     * A namespace for internal implementation details of the TrilinosWrapper
     * members.
     *
     * @ingroup TpetraWrappers
     */
    namespace internal
    {
      /**
       * Declare type for container size.
       */
      using size_type = types::global_dof_index;

      /**
       * This class implements a wrapper for accessing the Trilinos Tpetra
       * vector in the same way as we access deal.II objects: it is initialized
       * with a vector and an element within it, and has a conversion operator
       * to extract the scalar value of this element. It also has a variety of
       * assignment operator for writing to this one element.
       *
       * @ingroup TpetraWrappers
       */
      template <typename Number,
                typename MemorySpace = dealii::MemorySpace::Host>
      class VectorReference
      {
      private:
        /**
         * Constructor. It is made private so as to only allow the actual vector
         * class to create it.
         */
        VectorReference(Vector<Number, MemorySpace> &vector,
                        const size_type              index);

      public:
        /**
         * Copy constructor.
         */
        VectorReference(const VectorReference &) = default;

        /**
         * This looks like a copy operator, but does something different than
         * usual. In particular, it does not copy the member variables of this
         * reference. Rather, it handles the situation where we have two vectors
         * @p v and @p w, and assign elements like in <tt>v(i)=w(i)</tt>. Here,
         * both left and right hand side of the assignment have data type
         * VectorReference, but what we really mean is to assign the vector
         * elements represented by the two references. This operator implements
         * this operation. Note also that this allows us to make the assignment
         * operator const.
         */
        const VectorReference &
        operator=(const VectorReference &r) const;

        /**
         * Same as above but for non-const reference objects.
         */
        VectorReference &
        operator=(const VectorReference &r);

        /**
         * Set the referenced element of the vector to <tt>s</tt>.
         */
        const VectorReference &
        operator=(const Number &s) const;

        /**
         * Add <tt>s</tt> to the referenced element of the vector->
         */
        const VectorReference &
        operator+=(const Number &s) const;

        /**
         * Subtract <tt>s</tt> from the referenced element of the vector->
         */
        const VectorReference &
        operator-=(const Number &s) const;

        /**
         * Multiply the referenced element of the vector by <tt>s</tt>.
         */
        const VectorReference &
        operator*=(const Number &s) const;

        /**
         * Divide the referenced element of the vector by <tt>s</tt>.
         */
        const VectorReference &
        operator/=(const Number &s) const;

        /**
         * Convert the reference to an actual value, i.e. return the value of
         * the referenced element of the vector.
         */
        operator Number() const;

        /**
         * Exception
         */
        DeclException1(ExcTrilinosError,
                       int,
                       << "An error with error number " << arg1
                       << " occurred while calling a Trilinos function");

        /*
         * Access to a an element that is not (locally-)owned.
         *
         * @ingroup Exceptions
         */
        DeclException4(
          ExcAccessToNonLocalElement,
          size_type,
          size_type,
          size_type,
          size_type,
          << "You are trying to access element " << arg1
          << " of a distributed vector, but this element is not stored "
          << "on the current processor. Note: There are " << arg2
          << " elements stored "
          << "on the current processor from within the range [" << arg3 << ','
          << arg4 << "] but Trilinos vectors need not store contiguous "
          << "ranges on each processor, and not every element in "
          << "this range may in fact be stored locally."
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
         * Point to the vector we are referencing.
         */
        Vector<Number, MemorySpace> &vector;

        /**
         * Index of the referenced element of the vector.
         */
        const size_type index;

        // Make the vector class a friend, so that it can create objects of the
        // present type.
        friend class Vector<Number, MemorySpace>;
      }; // class VectorReference

    } // namespace internal
    /**
     * @endcond
     */


    /**
     * This class implements a wrapper to the Trilinos distributed vector
     * class Tpetra::Vector. This class requires Trilinos to be
     * compiled with MPI support.
     *
     * Moreover, this class takes an optional template argument for
     * the memory space used. By default, all memory is allocated on the CPU.
     *
     * @ingroup TpetraWrappers
     * @ingroup Vectors
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class Vector : public ReadVector<Number>
    {
    public:
      /**
       * Declare some of the standard types used in all containers.
       */
      using value_type = Number;
      using real_type  = typename numbers::NumberTraits<Number>::real_type;
      using size_type  = types::global_dof_index;
      using reference  = internal::VectorReference<Number, MemorySpace>;
      using const_reference =
        const internal::VectorReference<Number, MemorySpace>;

      /**
       * @name 1: Basic Object-handling
       */
      /** @{ */
      /**
       * Default constructor that generates an empty (zero size) vector. The
       * function <tt>reinit()</tt> will have to give the vector the correct
       * size and distribution among processes in case of an MPI run.
       */
      Vector();

      /**
       * Copy constructor. Sets the dimension and the partitioning to that of
       * the given vector and copies all elements.
       */
      Vector(const Vector &V);

      /**
       *  Copy constructor from Teuchos::RCP<Tpetra::Vector>.
       */
      Vector(
        const Teuchos::RCP<TpetraTypes::VectorType<Number, MemorySpace>> V);

      /**
       * TODO: This is not used
       * This constructor takes an IndexSet that defines how to distribute the
       * individual components among the MPI processors. Since it also
       * includes information about the size of the vector, this is all we
       * need to generate a %parallel vector.
       */
      explicit Vector(const IndexSet &parallel_partitioner,
                      const MPI_Comm  communicator);

      /**
       * In addition to just specifying one index set as in all the other
       * methods above, this method allows to supply an additional set of
       * ghost entries.
       *
       * Depending on whether the @p locally_relevant_or_ghost_entries argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"

       */
      explicit Vector(const IndexSet &locally_owned_entries,
                      const IndexSet &ghost_entries,
                      const MPI_Comm  communicator,
                      const bool      vector_writable = false);

      /**
       * Release all memory and return to a state just like after having called
       * the default constructor.
       */
      void
      clear();

      /**
       * Reinit functionality. This function destroys the old vector content
       * and generates a new one based on the input partitioning. The flag
       * <tt>omit_zeroing_entries</tt> determines whether the vector should be
       * filled with zeros (`false`) or left in an undetermined state (`true`).
       */
      void
      reinit(const IndexSet &parallel_partitioner,
             const MPI_Comm  communicator         = MPI_COMM_WORLD,
             const bool      omit_zeroing_entries = false);

      /**
       * Reinit functionality. This function destroys the old vector content
       * and generates a new one based on the input partitioning. In addition
       * to just specifying one index set as in all the other methods above,
       * this method allows to supply an additional set of ghost entries.
       *
       * Depending on whether the @p locally_relevant_or_ghost_entries argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      void
      reinit(const IndexSet &locally_owned_entries,
             const IndexSet &locally_relevant_or_ghost_entries,
             const MPI_Comm  communicator    = MPI_COMM_WORLD,
             const bool      vector_writable = false);

      /**
       * Change the dimension to that of the vector @p V. The elements of @p V are not
       * copied.
       */
      void
      reinit(const Vector<Number, MemorySpace> &V,
             const bool                         omit_zeroing_entries = false);

      /**
       * Swap the contents of this vector and the other vector @p v. One could do
       * this operation with a temporary variable and copying over the data
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
      virtual void
      swap(Vector &v) noexcept;

      /**
       * Extract a range of elements all at once.
       */
      virtual void
      extract_subvector_to(
        const ArrayView<const types::global_dof_index> &indices,
        const ArrayView<Number> &elements) const override;

      /**
       * Copy function. This function takes a Vector and copies all the
       * elements. The Vector will have the same parallel distribution as @p
       * V.
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
      Vector &
      operator=(const Vector &V);

      /**
       * Copy function. This function takes a Vector and copies all the
       * elements. The Vector will have the same parallel distribution as @p
       * V.
       */
      template <typename OtherNumber>
      Vector &
      operator=(const dealii::Vector<OtherNumber> &V);

      /**
       * Sets all elements of the vector to the scalar @p s. This operation is
       * only allowed if @p s is equal to zero.
       */
      Vector &
      operator=(const Number s);

      /**
       * Imports all the elements present in the vector's IndexSet from the
       * input vector @p V. VectorOperation::values @p operation is used to decide if
       * the elements in @p V should be added to the current vector or replace the
       * current elements. The last parameter can be used if the same
       * communication pattern is used multiple times. This can be used to
       * improve performance.
       */
      void
      import_elements(
        const ReadWriteVector<Number> &V,
        VectorOperation::values        operation,
        const Teuchos::RCP<const Utilities::MPI::CommunicationPatternBase>
          &communication_pattern);

      /**
       * @deprecated Use Teuchos::RCP<> instead of std::shared_ptr<>.
       */
      DEAL_II_DEPRECATED
      void
      import_elements(
        const ReadWriteVector<Number> &V,
        VectorOperation::values        operation,
        const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
          &communication_pattern);

      /*
       * Imports all the elements present in the vector's IndexSet from the
       * input vector @p V. VectorOperation::values @p operation is used to decide if
       * the elements in @p V should be added to the current vector or replace the
       * current elements.
       */
      void
      import_elements(const ReadWriteVector<Number> &V,
                      VectorOperation::values        operation);

      /**
       * @deprecated Use import_elements() instead.
       */
      DEAL_II_DEPRECATED
      void
      import(const ReadWriteVector<Number> &V,
             VectorOperation::values        operation,
             std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
               communication_pattern = {})
      {
        import_elements(V, operation, communication_pattern);
      }

      /** @} */


      /**
       * @name 2: Data-Access
       */
      /** @{ */

      /**
       * Provide access to a given element, both read and write.
       *
       * When using a vector distributed with MPI, this operation only makes
       * sense for elements that are actually present on the calling processor.
       * Otherwise, an exception is thrown.
       */
      reference
      operator()(const size_type index);

      /**
       * Provide read-only access to an element.
       *
       * When using a vector distributed with MPI, this operation only makes
       * sense for elements that are actually present on the calling processor.
       * Otherwise, an exception is thrown.
       */
      Number
      operator()(const size_type index) const;

      /**
       * Provide access to a given element, both read and write.
       *
       * Exactly the same as operator().
       */
      reference
      operator[](const size_type index);

      /**
       * Provide read-only access to an element.
       *
       * Exactly the same as operator().
       */
      Number
      operator[](const size_type index) const;

      /** @} */


      /**
       * @name 3: Modification of vectors
       */
      /** @{ */

      /**
       * Multiply the entire vector by a fixed factor.
       */
      Vector &
      operator*=(const Number factor);

      /**
       * Divide the entire vector by a fixed factor.
       */
      Vector &
      operator/=(const Number factor);

      /**
       * Add the vector @p V to the present one.
       */
      Vector &
      operator+=(const Vector<Number, MemorySpace> &V);

      /**
       * Subtract the vector @p V from the present one.
       */
      Vector &
      operator-=(const Vector<Number, MemorySpace> &V);

      /**
       * Return the scalar product of two vectors. The vectors need to have the
       * same layout.
       */
      Number
      operator*(const Vector<Number, MemorySpace> &V) const;

      /**
       * Add @p a to all components. Note that @p is a scalar not a vector.
       */
      void
      add(const Number a);

      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this +=
       * a*V</tt>. The vectors need to have the same layout.
       */
      void
      add(const Number a, const Vector<Number, MemorySpace> &V);

      /**
       * Multiple addition of multiple of a vector, i.e. <tt>*this> +=
       * a*V+b*W</tt>. The vectors need to have the same layout.
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
       * Take an address where <tt>n_elements</tt> are stored contiguously and
       * add them into the vector. Handles all cases which are not covered by
       * the other two <tt>add()</tt> functions above.
       */
      void
      add(const size_type  n_elements,
          const size_type *indices,
          const Number    *values);

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this
       * = s*(*this)+a*V</tt>.
       */
      void
      sadd(const Number                       s,
           const Number                       a,
           const Vector<Number, MemorySpace> &V);

      /**
       * A collective set operation: instead of setting individual elements of a
       * vector, this function allows to set a whole set of elements at once.
       * It is assumed that the elements to be set are located in contiguous
       * memory.
       */
      void
      set(const size_type  n_elements,
          const size_type *indices,
          const Number    *values);

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument. This function is mostly meant to simulate multiplication
       * (and immediate re-assignment) by a diagonal scaling matrix. The
       * vectors need to have the same layout.
       */
      void
      scale(const Vector<Number, MemorySpace> &scaling_factors);

      /**
       * Assignment <tt>*this = a*V</tt>.
       */
      void
      equ(const Number a, const Vector<Number, MemorySpace> &V);

      /**
       * Return whether the vector contains only elements with value zero.
       */
      bool
      all_zero() const;

      /**
       * Return @p true if the vector has no negative entries, i.e. all entries
       * are zero or positive. This function is used, for example, to check
       * whether refinement indicators are really all positive (or zero).
       */
      bool
      is_non_negative() const;

      /** @} */


      /**
       * @name 4: Scalar products, norms and related operations
       */
      /** @{ */

      /**
       * Return the mean value of the element of this vector.
       */
      Number
      mean_value() const;

      /**
       * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
       * absolute values of all entries among all processors).
       */
      real_type
      l1_norm() const;

      /**
       * Return the l<sub>2</sub> norm of the vector (i.e., the square root of
       * the sum of the square of all entries among all processors).
       */
      real_type
      l2_norm() const;

      /**
       * Return the maximum norm of the vector (i.e., the maximum absolute value
       * among all entries and among all processors).
       */
      real_type
      linfty_norm() const;

      /**
       * Return the square of the l<sub>2</sub>-norm.
       */
      real_type
      norm_sqr() const;

      /**
       *
       *
       */
      /**
       * Performs a combined operation of a vector addition and a subsequent
       * inner product, returning the value of the inner product. In other
       * words, the result of this function is the same as if the user called
       * @code
       * this->add(a, V);
       * return_value = *this * W;
       * @endcode
       *
       * The reason this function exists is that this operation involves less
       * memory transfer than calling the two functions separately. This
       * method only needs to load three vectors, @p this, @p V, @p W, whereas
       * calling separate methods means to load the calling vector @p this
       * twice. Since most vector operations are memory transfer limited, this
       * reduces the time by 25\% (or 50\% if @p W equals @p this).
       *
       * The vectors need to have the same layout.
       *
       * For complex-valued vectors, the scalar product in the second step is
       * implemented as
       * $\left<v,w\right>=\sum_i v_i \bar{w_i}$.
       */
      Number
      add_and_dot(const Number                       a,
                  const Vector<Number, MemorySpace> &V,
                  const Vector<Number, MemorySpace> &W);

      /** @} */


      /**
       * @name 5: Scalar products, norms and related operations
       */
      /** @{ */

      /**
       * Return whether the vector has ghost elements or not.
       */
      bool
      has_ghost_elements() const;

      /**
       * Test for equality. This function assumes that the present vector and
       * the one to compare with have the same size already, since comparing
       * vectors of different sizes makes not much sense anyway.
       */
      bool
      operator==(const Vector<Number, MemorySpace> &v) const;

      /**
       * Test for inequality. This function assumes that the present vector and
       * the one to compare with have the same size already, since comparing
       * vectors of different sizes makes not much sense anyway.
       */
      bool
      operator!=(const Vector<Number, MemorySpace> &v) const;

      /**
       * Return the global size of the vector, equal to the sum of the number of
       * locally owned indices among all processors.
       */
      virtual size_type
      size() const override;

      /**
       * Return the local size of the vector, i.e., the number of indices
       * owned locally.
       */
      size_type
      locally_owned_size() const;

      /**
       * Return a pair of indices indicating which elements of this vector are
       * stored locally. The first number is the index of the first element
       * stored, the second the index of the one past the last one that is
       * stored locally. If this is a sequential vector, then the result will be
       * the pair <code>(0,N)</code>, otherwise it will be a pair
       * <code>(i,i+n)</code>, where <code>n</code> is the number of elements
       * stored on this processor and and <code>i</code> is the first element of
       * the vector stored on this processor, corresponding to the half open
       * interval $[i,i+n)$
       *
       * @note The description above is true most of the time, but not always.
       * In particular, Trilinos vectors need not store contiguous ranges of
       * elements such as $[i,i+n)$. Rather, it can store vectors where the
       * elements are distributed in an arbitrary way across all processors and
       * each processor simply stores a particular subset, not necessarily
       * contiguous. In this case, this function clearly makes no sense since it
       * could, at best, return a range that includes all elements that are
       * stored locally. Thus, the function only succeeds if the locally stored
       * range is indeed contiguous. It will trigger an assertion if the local
       * portion of the vector is not contiguous.
       */
      std::pair<size_type, size_type>
      local_range() const;

      /**
       * Return whether @p index is in the local range or not, see also
       * local_range().
       *
       * @note The same limitation for the applicability of this function
       * applies as listed in the documentation of local_range().
       */
      bool
      in_local_range(const size_type index) const;

      /**
       * Return the state of the vector, i.e., whether compress() needs to be
       * called after an operation requiring data exchange. A call to compress()
       * is also needed when the method set() or add() has been called.
       */
      bool
      is_compressed() const;

      /**
       * Return the underlying MPI communicator.
       */
      MPI_Comm
      get_mpi_communicator() const;

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
      ::dealii::IndexSet
      locally_owned_elements() const;

      /** @} */


      /**
       * @name 6: Mixed stuff
       */
      /** @{ */
      /**
       * Compress the underlying representation of the Trilinos object, i.e.
       * flush the buffers of the vector object if it has any. This function is
       * necessary after writing into a vector element-by-element and before
       * anything else can be done on it.
       *
       * @param operation The compress mode (<code>Add</code> or <code>Insert</code>)
       * in case the vector has not been written to since the last time this
       * function was called. The argument is ignored if the vector has been
       * added or written to since the last time compress() was called.
       *
       * See
       * @ref GlossCompress "Compressing distributed objects"
       * for more information.
       */
      void
      compress(const VectorOperation::values operation);

      /**
       * Return a const reference to the underlying Trilinos
       * Tpetra::Vector class.
       */
      const TpetraTypes::VectorType<Number, MemorySpace> &
      trilinos_vector() const;

      /**
       * Return a (modifiable) reference to the underlying Trilinos
       * Tpetra::Vector class.
       */
      TpetraTypes::VectorType<Number, MemorySpace> &
      trilinos_vector();

      /**
       * Return a const Teuchos::RCP to the underlying Trilinos
       * Tpetra::Vector class.
       */
      Teuchos::RCP<const TpetraTypes::VectorType<Number, MemorySpace>>
      trilinos_rcp() const;

      /**
       * Return a (modifiable) Teuchos::RCP to the underlying Trilinos
       * Tpetra::Vector class.
       */
      Teuchos::RCP<TpetraTypes::VectorType<Number, MemorySpace>>
      trilinos_rcp();

      /**
       * Prints the vector to the output stream @p out.
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

      /**
       * Return the mpi communicator
       */
      MPI_Comm
      mpi_comm() const;

      /** @} */


      /**
       * The vectors have different partitioning, i.e. their IndexSet objects
       * don't represent the same indices.
       *
       * @ingroup Exceptions
       */
      DeclException0(ExcDifferentParallelPartitioning);

      /**
       * Attempt to perform an operation between two incompatible vector types.
       *
       * @ingroup Exceptions
       */
      DeclException0(ExcVectorTypeNotCompatible);

      /*
       * Access to a an element that is not (locally-)owned.
       *
       * @ingroup Exceptions
       */
      DeclException4(
        ExcAccessToNonLocalElement,
        size_type,
        size_type,
        size_type,
        size_type,
        << "You are trying to access element " << arg1
        << " of a distributed vector, but this element is not stored "
        << "on the current processor. Note: There are " << arg2
        << " elements stored "
        << "on the current processor from within the range [" << arg3 << ','
        << arg4 << "] but Trilinos vectors need not store contiguous "
        << "ranges on each processor, and not every element in "
        << "this range may in fact be stored locally."
        << "\n\n"
        << "A common source for this kind of problem is that you "
        << "are passing a 'fully distributed' vector into a function "
        << "that needs read access to vector elements that correspond "
        << "to degrees of freedom on ghost cells (or at least to "
        << "'locally active' degrees of freedom that are not also "
        << "'locally owned'). You need to pass a vector that has these "
        << "elements as ghost entries.");

      /**
       * Missing index set.
       *
       * @ingroup Exceptions
       */
      DeclExceptionMsg(ExcMissingIndexSet,
                       "To compress a vector, a locally_relevant_dofs "
                       "index set, and a locally_owned_dofs index set "
                       "must be provided. These index sets must be "
                       "provided either when the vector is initialized "
                       "or when compress is called. See the documentation "
                       "of compress() for more information.");

      /**
       * Exception thrown by an error in Trilinos.
       *
       * @ingroup Exceptions
       */
      DeclException1(ExcTrilinosError,
                     int,
                     << "An error with error number " << arg1
                     << " occurred while calling a Trilinos function");

    private:
      /**
       * Create the CommunicationPattern for the communication between the
       * IndexSet @p source_index_set and the current vector based
       * on the communicator @p mpi_comm.
       */
      void
      create_tpetra_comm_pattern(const IndexSet &source_index_set,
                                 const MPI_Comm  mpi_comm);

      /**
       * A boolean variable to hold information on whether the vector is
       * compressed or not.
       */
      bool compressed;

      /**
       * Store whether the vector has ghost elements or not.
       *
       * If the vector has no ghost elements, it can only access and modify
       * entries included in the locally owned index set.
       * And if the vector has ghost elements it can access and modify
       * entries included in the locally relevant index set.
       */
      bool has_ghost;

      /**
       * Teuchos::RCP to the actual Tpetra vector object.
       */
      Teuchos::RCP<TpetraTypes::VectorType<Number, MemorySpace>> vector;

      /**
       * A vector object in Trilinos to be used for collecting the non-local
       * elements if the vector was constructed with an additional IndexSet
       * describing ghost elements.
       */
      Teuchos::RCP<TpetraTypes::VectorType<Number, MemorySpace>>
        nonlocal_vector;

      /**
       * IndexSet of the elements of the last imported vector.
       */
      dealii::IndexSet source_stored_elements;

      /**
       * CommunicationPattern for the communication between the
       * source_stored_elements IndexSet and the current vector.
       */
      Teuchos::RCP<const TpetraWrappers::CommunicationPattern<MemorySpace>>
        tpetra_comm_pattern;

      // Make the reference class a friend.
      friend class internal::VectorReference<Number, MemorySpace>;
    };


    /* ------------------------- Inline functions ---------------------- */

    template <typename Number, typename MemorySpace>
    inline void
    swap(Vector<Number, MemorySpace> &u,
         Vector<Number, MemorySpace> &v) noexcept
    {
      u.swap(v);
    }


    template <typename Number, typename MemorySpace>
    inline bool
    Vector<Number, MemorySpace>::has_ghost_elements() const
    {
      return has_ghost;
    }



    template <typename Number, typename MemorySpace>
    inline bool
    Vector<Number, MemorySpace>::is_compressed() const
    {
      return compressed;
    }

    template <typename Number, typename MemorySpace>
    inline void
    Vector<Number, MemorySpace>::swap(Vector<Number, MemorySpace> &v) noexcept
    {
      vector.swap(v.vector);
    }


    template <typename Number, typename MemorySpace>
    inline void
    Vector<Number, MemorySpace>::add(const std::vector<size_type> &indices,
                                     const std::vector<Number>    &values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      AssertDimension(indices.size(), values.size());

      add(indices.size(), indices.data(), values.data());
    }



    template <typename Number, typename MemorySpace>
    inline void
    Vector<Number, MemorySpace>::add(const size_type  n_elements,
                                     const size_type *indices,
                                     const Number    *values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      auto vector_2d_local = vector->template getLocalView<Kokkos::HostSpace>(
        Tpetra::Access::ReadWrite);

      // Having extracted a view into the multivectors above, now also
      // extract a view into the one vector we actually store. We can
      // do this right away for the locally owned part. We defer creating
      // the view into the nonlocal part to when we know that we actually
      // need it; this also makes sure that we correctly deal with the
      // case where we do not actually store a nonlocal part.
      auto vector_1d_local = Kokkos::subview(vector_2d_local, Kokkos::ALL(), 0);
      using ViewType1d     = decltype(vector_1d_local);
      std::optional<ViewType1d> vector_1d_nonlocal;

#  if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      // Mark vector as to-be-modified. We may do the same with
      // the nonlocal part too if we end up writing into it.
      vector->template modify<Kokkos::HostSpace>();
#  endif

      for (size_type i = 0; i < n_elements; ++i)
        {
          const size_type row = indices[i];

          // Check if the index is in the locally owned index set.
          // If so, we can write right into the locally owned
          // part of the vector.
          if (TrilinosWrappers::types::int_type local_row =
                vector->getMap()->getLocalElement(row);
              local_row != Teuchos::OrdinalTraits<int>::invalid())
            {
              vector_1d_local(local_row) += values[i];

              // Set the compressed state to false only if there is nonlocal
              // part in this distributed vector, otherwise it's always
              // compressed.
              if (nonlocal_vector.get() != nullptr)
                compressed = false;
            }
          else
            {
              // If the element was not in the locally owned part,
              // we need to figure out whether it is in the nonlocal
              // part. It better be:
              Assert(nonlocal_vector.get() != nullptr, ExcInternalError());
              TrilinosWrappers::types::int_type nonlocal_row =
                nonlocal_vector->getMap()->getLocalElement(row);

#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
              Assert(nonlocal_row != Teuchos::OrdinalTraits<int>::invalid(),
                     ExcAccessToNonLocalElement(
                       row,
                       vector->getMap()->getLocalNumElements(),
                       vector->getMap()->getMinLocalIndex(),
                       vector->getMap()->getMaxLocalIndex()));
#  else
              Assert(nonlocal_row != Teuchos::OrdinalTraits<int>::invalid(),
                     ExcAccessToNonLocalElement(
                       row,
                       vector->getMap()->getNodeNumElements(),
                       vector->getMap()->getMinLocalIndex(),
                       vector->getMap()->getMaxLocalIndex()));

#  endif

              // Having asserted that it is, write into the nonlocal part.
              // To do so, we first need to make sure that we have a view
              // of the nonlocal part of the vectors, since we have
              // deferred creating this view previously:
              if (!vector_1d_nonlocal)
                {
#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
                  auto vector_2d_nonlocal =
                    nonlocal_vector->template getLocalView<Kokkos::HostSpace>(
                      Tpetra::Access::ReadWrite);
#  else
                  auto vector_2d_nonlocal =
                    nonlocal_vector->template getLocalView<Kokkos::HostSpace>();
#  endif

                  vector_1d_nonlocal =
                    Kokkos::subview(vector_2d_nonlocal, Kokkos::ALL(), 0);

#  if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
                  // Mark the nonlocal vector as to-be-modified as well.
                  nonlocal_vector->template modify<Kokkos::HostSpace>();
#  endif
                }
              (*vector_1d_nonlocal)(nonlocal_row) += values[i];
              compressed = false;
            }
        }

#  if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      vector->template sync<
        typename Tpetra::Vector<Number, int, types::signed_global_dof_index>::
          device_type::memory_space>();

      // If we have created a view to the nonlocal part, then we have also
      // written into it. Flush these modifications.
      if (vector_1d_nonlocal)
        nonlocal_vector->template sync<
          typename Tpetra::Vector<Number, int, types::signed_global_dof_index>::
            device_type::memory_space>();
#  endif
    }



    template <typename Number, typename MemorySpace>
    inline void
    Vector<Number, MemorySpace>::set(const size_type  n_elements,
                                     const size_type *indices,
                                     const Number    *values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      auto vector_2d_local = vector->template getLocalView<Kokkos::HostSpace>(
        Tpetra::Access::ReadWrite);
#  else
      vector->template sync<Kokkos::HostSpace>();

      auto vector_2d_local = vector->template getLocalView<Kokkos::HostSpace>();
#  endif

      // Having extracted a view into the multivectors above, now also
      // extract a view into the one vector we actually store. We can
      // do this right away for the locally owned part. We defer creating
      // the view into the nonlocal part to when we know that we actually
      // need it; this also makes sure that we correctly deal with the
      // case where we do not actually store a nonlocal part.
      auto vector_1d_local = Kokkos::subview(vector_2d_local, Kokkos::ALL(), 0);
      using ViewType1d     = decltype(vector_1d_local);
      std::optional<ViewType1d> vector_1d_nonlocal;

#  if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      // Mark vector as to-be-modified. We may do the same with
      // the nonlocal part too if we end up writing into it.
      vector->template modify<Kokkos::HostSpace>();
#  endif

      for (size_type i = 0; i < n_elements; ++i)
        {
          const size_type row = indices[i];

          // Check if the index is in the locally owned index set.
          // If so, we can write right into the locally owned
          // part of the vector.
          if (TrilinosWrappers::types::int_type local_row =
                vector->getMap()->getLocalElement(row);
              local_row != Teuchos::OrdinalTraits<int>::invalid())
            {
              vector_1d_local(local_row) = values[i];

              // Set the compressed state to false only if there is nonlocal
              // part in this distributed vector, otherwise it's always
              // compressed.
              if (nonlocal_vector.get() != nullptr)
                compressed = false;
            }
          else
            {
              // If the element was not in the locally owned part,
              // we need to figure out whether it is in the nonlocal
              // part. It better be:
              Assert(nonlocal_vector.get() != nullptr, ExcInternalError());
              TrilinosWrappers::types::int_type nonlocal_row =
                nonlocal_vector->getMap()->getLocalElement(row);

#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
              Assert(nonlocal_row != Teuchos::OrdinalTraits<int>::invalid(),
                     ExcAccessToNonLocalElement(
                       row,
                       vector->getMap()->getLocalNumElements(),
                       vector->getMap()->getMinLocalIndex(),
                       vector->getMap()->getMaxLocalIndex()));
#  else
              Assert(nonlocal_row != Teuchos::OrdinalTraits<int>::invalid(),
                     ExcAccessToNonLocalElement(
                       row,
                       vector->getMap()->getNodeNumElements(),
                       vector->getMap()->getMinLocalIndex(),
                       vector->getMap()->getMaxLocalIndex()));

#  endif

              // Having asserted that it is, write into the nonlocal part.
              // To do so, we first need to make sure that we have a view
              // of the nonlocal part of the vectors, since we have
              // deferred creating this view previously:
              if (!vector_1d_nonlocal)
                {
#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
                  auto vector_2d_nonlocal =
                    nonlocal_vector->template getLocalView<Kokkos::HostSpace>(
                      Tpetra::Access::ReadWrite);
#  else
                  auto vector_2d_nonlocal =
                    nonlocal_vector->template getLocalView<Kokkos::HostSpace>();
#  endif

                  vector_1d_nonlocal =
                    Kokkos::subview(vector_2d_nonlocal, Kokkos::ALL(), 0);

#  if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
                  // Mark the nonlocal vector as to-be-modified as well.
                  nonlocal_vector->template modify<Kokkos::HostSpace>();
#  endif
                }
              (*vector_1d_nonlocal)(nonlocal_row) = values[i];
              compressed                          = false;
            }
        }

#  if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      vector->template sync<
        typename Tpetra::Vector<Number, int, types::signed_global_dof_index>::
          device_type::memory_space>();

      // If we have created a view to the nonlocal part, then we have also
      // written into it. Flush these modifications.
      if (vector_1d_nonlocal)
        nonlocal_vector->template sync<
          typename Tpetra::Vector<Number, int, types::signed_global_dof_index>::
            device_type::memory_space>();
#  endif
    }



    template <typename Number, typename MemorySpace>
    inline internal::VectorReference<Number, MemorySpace>
    Vector<Number, MemorySpace>::operator()(const size_type index)
    {
      return internal::VectorReference(*this, index);
    }

    template <typename Number, typename MemorySpace>
    inline internal::VectorReference<Number, MemorySpace>
    Vector<Number, MemorySpace>::operator[](const size_type index)
    {
      return operator()(index);
    }

    template <typename Number, typename MemorySpace>
    inline Number
    Vector<Number, MemorySpace>::operator[](const size_type index) const
    {
      return operator()(index);
    }

#  ifndef DOXYGEN

    // VectorReference
    namespace internal
    {
      template <typename Number, typename MemorySpace>
      inline VectorReference<Number, MemorySpace>::VectorReference(
        Vector<Number, MemorySpace> &vector,
        const size_type              index)
        : vector(vector)
        , index(index)
      {}



      template <typename Number, typename MemorySpace>
      inline const VectorReference<Number, MemorySpace> &
      VectorReference<Number, MemorySpace>::operator=(
        const VectorReference<Number, MemorySpace> &r) const
      {
        // as explained in the class
        // documentation, this is not the copy
        // operator. so simply pass on to the
        // "correct" assignment operator
        *this = static_cast<Number>(r);

        return *this;
      }



      template <typename Number, typename MemorySpace>
      inline VectorReference<Number, MemorySpace> &
      VectorReference<Number, MemorySpace>::operator=(
        const VectorReference<Number, MemorySpace> &r)
      {
        // as above
        *this = static_cast<Number>(r);

        return *this;
      }



      template <typename Number, typename MemorySpace>
      inline const VectorReference<Number, MemorySpace> &
      VectorReference<Number, MemorySpace>::operator=(const Number &value) const
      {
        vector.set(1, &index, &value);
        return *this;
      }



      template <typename Number, typename MemorySpace>
      inline const VectorReference<Number, MemorySpace> &
      VectorReference<Number, MemorySpace>::operator+=(
        const Number &value) const
      {
        vector.add(1, &index, &value);
        return *this;
      }



      template <typename Number, typename MemorySpace>
      inline const VectorReference<Number, MemorySpace> &
      VectorReference<Number, MemorySpace>::operator-=(
        const Number &value) const
      {
        Number new_value = -value;
        vector.add(1, &index, &new_value);
        return *this;
      }



      template <typename Number, typename MemorySpace>
      inline const VectorReference<Number, MemorySpace> &
      VectorReference<Number, MemorySpace>::operator*=(
        const Number &value) const
      {
        Number new_value = static_cast<Number>(*this) * value;
        vector.set(1, &index, &new_value);
        return *this;
      }



      template <typename Number, typename MemorySpace>
      inline const VectorReference<Number, MemorySpace> &
      VectorReference<Number, MemorySpace>::operator/=(
        const Number &value) const
      {
        Number new_value = static_cast<Number>(*this) / value;
        vector.set(1, &index, &new_value);
        return *this;
      }
    } // namespace internal

#  endif /* DOXYGEN */

  } // namespace TpetraWrappers

  /** @} */

} // namespace LinearAlgebra


namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename>
    class ReinitHelper;

    /**
     * A helper class used internally in linear_operator.h. Specialization for
     * LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>.
     */
    template <typename Number, typename MemorySpace>
    class ReinitHelper<
      LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>>
    {
    public:
      template <typename Matrix>
      static void
      reinit_range_vector(
        const Matrix                                               &matrix,
        LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace> &v,
        bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_range_indices(),
                 matrix.get_mpi_communicator(),
                 omit_zeroing_entries);
      }

      template <typename Matrix>
      static void
      reinit_domain_vector(
        const Matrix                                               &matrix,
        LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace> &v,
        bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_domain_indices(),
                 matrix.get_mpi_communicator(),
                 omit_zeroing_entries);
      }
    };

  } // namespace LinearOperatorImplementation
} /* namespace internal */

/**
 * Declare dealii::LinearAlgebra::TpetraWrappers::Vector as distributed vector.
 */
template <typename Number, typename MemorySpace>
struct is_serial_vector<
  LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace>> : std::false_type
{};

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
