// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_epetra_vector_h
#define dealii_trilinos_epetra_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/mpi_stub.h>

#  include <deal.II/lac/read_vector.h>
#  include <deal.II/lac/trilinos_epetra_communication_pattern.h>
#  include <deal.II/lac/vector_operation.h>
#  include <deal.II/lac/vector_type_traits.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <Epetra_FEVector.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  // Forward declaration
  template <typename Number>
  class ReadWriteVector;

  /**
   * A namespace for classes that provide wrappers for Trilinos' Epetra vectors.
   *
   * This namespace provides wrappers for the Epetra_FEVector class from the
   * Epetra package (https://trilinos.github.io/epetra.html) that is part of
   * Trilinos.
   */
  namespace EpetraWrappers
  {
    /**
     * This class defines type aliases that are used in vector classes
     * within the EpetraWrappers namespace.
     */
    class VectorTraits
    {
    public:
      using value_type = double;
      using size_type  = dealii::types::global_dof_index;
    };

    /**
     * @cond internal
     */

    /**
     * A namespace for internal implementation details of the TrilinosWrapper
     * members.
     *
     * @ingroup TrilinosWrappers
     */
    namespace internal
    {
      /**
       * This class implements a wrapper for accessing the Trilinos Epetra
       * vector in the same way as we access deal.II objects: it is initialized
       * with a vector and an element within it, and has a conversion operator
       * to extract the scalar value of this element. It also has a variety of
       * assignment operator for writing to this one element.
       *
       * @ingroup TrilinosWrappers
       */
      class VectorReference
      {
      private:
        using value_type = VectorTraits::value_type;
        using size_type  = VectorTraits::size_type;

        /**
         * Constructor. It is made private so as to only allow the actual vector
         * class to create it.
         */
        VectorReference(Vector &vector, const size_type index);

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
        operator=(const value_type &s) const;

        /**
         * Add <tt>s</tt> to the referenced element of the vector->
         */
        const VectorReference &
        operator+=(const value_type &s) const;

        /**
         * Subtract <tt>s</tt> from the referenced element of the vector->
         */
        const VectorReference &
        operator-=(const value_type &s) const;

        /**
         * Multiply the referenced element of the vector by <tt>s</tt>.
         */
        const VectorReference &
        operator*=(const value_type &s) const;

        /**
         * Divide the referenced element of the vector by <tt>s</tt>.
         */
        const VectorReference &
        operator/=(const value_type &s) const;

        /**
         * Convert the reference to an actual value, i.e. return the value of
         * the referenced element of the vector.
         */
        operator value_type() const;

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
        Vector &vector;

        /**
         * Index of the referenced element of the vector.
         */
        const size_type index;

        // Make the vector class a friend, so that it can create objects of the
        // present type.
        friend class ::dealii::LinearAlgebra::EpetraWrappers::Vector;
      }; // class VectorReference

    } // namespace internal
    /**
     * @endcond
     */


    /**
     * This class implements a wrapper to the Trilinos distributed vector
     * class Epetra_FEVector. This class requires Trilinos to be compiled
     * with MPI support.
     *
     * @ingroup TrilinosWrappers
     * @ingroup Vectors
     */
    class Vector : public ReadVector<VectorTraits::value_type>
    {
    public:
      using value_type      = VectorTraits::value_type;
      using size_type       = VectorTraits::size_type;
      using real_type       = value_type;
      using reference       = internal::VectorReference;
      using const_reference = const internal::VectorReference;

      /**
       * @name 1: Basic Object-handling
       */
      /** @{ */

      /**
       * Constructor. Create a vector of dimension zero.
       */
      Vector();

      /**
       * Copy constructor. Sets the dimension and the partitioning to that of
       * the given vector and copies all elements.
       */
      Vector(const Vector &V);

      /**
       * This constructor takes an IndexSet that defines how to distribute the
       * individual components among the MPI processors. Since it also
       * includes information about the size of the vector, this is all we
       * need to generate a %parallel vector.
       */
      explicit Vector(const IndexSet &parallel_partitioner,
                      const MPI_Comm  communicator);

      /**
       * Reinit functionality. This function destroys the old vector content
       * and generates a new one based on the input partitioning. The flag
       * <tt>omit_zeroing_entries</tt> determines whether the vector should be
       * filled with zeros (`false`) or left in an undetermined state (`true`).
       */
      void
      reinit(const IndexSet &parallel_partitioner,
             const MPI_Comm  communicator,
             const bool      omit_zeroing_entries = false);

      /**
       * Change the dimension to that of the vector @p V. The elements of @p V are not
       * copied.
       */
      void
      reinit(const Vector &V, const bool omit_zeroing_entries = false);

      /**
       * Extract a range of elements all at once.
       */
      virtual void
      extract_subvector_to(
        const ArrayView<const types::global_dof_index> &indices,
        const ArrayView<double> &elements) const override;

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
       * Sets all elements of the vector to the scalar @p s. This operation is
       * only allowed if @p s is equal to zero.
       */
      Vector &
      operator=(const double s);

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
        const ReadWriteVector<double> &V,
        VectorOperation::values        operation,
        const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
          &communication_pattern = {});

      /**
       * @deprecated Use import_elements() instead.
       */
      DEAL_II_DEPRECATED
      void
      import(
        const ReadWriteVector<double> &V,
        const VectorOperation::values  operation,
        const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
          &communication_pattern = {})
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
      value_type
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
      value_type
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
      operator*=(const double factor);

      /**
       * Divide the entire vector by a fixed factor.
       */
      Vector &
      operator/=(const double factor);

      /**
       * Add the vector @p V to the present one.
       */
      Vector &
      operator+=(const Vector &V);

      /**
       * Subtract the vector @p V from the present one.
       */
      Vector &
      operator-=(const Vector &V);

      /**
       * Return the scalar product of two vectors. The vectors need to have the
       * same layout.
       */
      double
      operator*(const Vector &V) const;

      /**
       * Add @p a to all components. Note that @p is a scalar not a vector.
       */
      void
      add(const double a);

      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this +=
       * a*V</tt>. The vectors need to have the same layout.
       */
      void
      add(const double a, const Vector &V);

      /**
       * Multiple addition of multiple of a vector, i.e. <tt>*this> +=
       * a*V+b*W</tt>. The vectors need to have the same layout.
       */
      void
      add(const double a, const Vector &V, const double b, const Vector &W);

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this
       * = s*(*this)+a*V</tt>.
       */
      void
      sadd(const double s, const double a, const Vector &V);

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument. This function is mostly meant to simulate multiplication
       * (and immediate re-assignment) by a diagonal scaling matrix. The
       * vectors need to have the same layout.
       */
      void
      scale(const Vector &scaling_factors);

      /**
       * Assignment <tt>*this = a*V</tt>.
       */
      void
      equ(const double a, const Vector &V);

      /**
       * Return whether the vector contains only elements with value zero.
       */
      bool
      all_zero() const;

      /** @} */


      /**
       * @name 4: Scalar products, norms and related operations
       */
      /** @{ */

      /**
       * Return the mean value of the element of this vector.
       */
      double
      mean_value() const;

      /**
       * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
       * absolute values of all entries among all processors).
       */
      double
      l1_norm() const;

      /**
       * Return the l<sub>2</sub> norm of the vector (i.e., the square root of
       * the sum of the square of all entries among all processors).
       */
      double
      l2_norm() const;

      /**
       * Return the maximum norm of the vector (i.e., the maximum absolute value
       * among all entries and among all processors).
       */
      double
      linfty_norm() const;

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
      double
      add_and_dot(const double a, const Vector &V, const Vector &W);

      /** @} */

      /**
       * @name 5: Scalar products, norms and related operations
       */
      /** @{ */

      /**
       * This function always returns false and is present only for backward
       * compatibility.
       */
      bool
      has_ghost_elements() const;

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
       * See
       * @ref GlossCompress "Compressing distributed objects"
       * for more information.
       */
      void
      compress(const VectorOperation::values operation);

      /**
       * Return a const reference to the underlying Trilinos
       * Epetra_FEVector class.
       */
      const Epetra_FEVector &
      trilinos_vector() const;

      /**
       * Return a (modifiable) reference to the underlying Trilinos
       * Epetra_FEVector class.
       */
      Epetra_FEVector &
      trilinos_vector();

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

      /** @} */

      /**
       * The vectors have different partitioning, i.e. their IndexSet objects
       * don't represent the same indices.
       */
      DeclException0(ExcDifferentParallelPartitioning);

      /**
       * Attempt to perform an operation between two incompatible vector types.
       *
       * @ingroup Exceptions
       */
      DeclException0(ExcVectorTypeNotCompatible);

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
      create_epetra_comm_pattern(const IndexSet &source_index_set,
                                 const MPI_Comm  mpi_comm);

      /**
       * Pointer to the actual Epetra vector object.
       */
      std::unique_ptr<Epetra_FEVector> vector;

      /**
       * IndexSet of the elements of the last imported vector.
       */
      ::dealii::IndexSet source_stored_elements;

      /**
       * CommunicationPattern for the communication between the
       * source_stored_elements IndexSet and the current vector.
       */
      std::shared_ptr<const CommunicationPattern> epetra_comm_pattern;

      // Make the reference class a friend.
      friend class internal::VectorReference;
    };

#  ifndef DOXYGEN

    // VectorReference
    namespace internal
    {
      inline VectorReference::VectorReference(Vector         &vector,
                                              const size_type index)
        : vector(vector)
        , index(index)
      {}



      inline const VectorReference &
      VectorReference::operator=(const VectorReference &r) const
      {
        // as explained in the class
        // documentation, this is not the copy
        // operator. so simply pass on to the
        // "correct" assignment operator
        *this = static_cast<value_type>(r);

        return *this;
      }



      inline VectorReference &
      VectorReference::operator=(const VectorReference &r)
      {
        // as above
        *this = static_cast<value_type>(r);

        return *this;
      }



      inline const VectorReference &
      VectorReference::operator=(const value_type &value) const
      {
        (*vector.vector)[0][index] = value;

        return *this;
      }



      inline const VectorReference &
      VectorReference::operator+=(const value_type &value) const
      {
        value_type new_value       = static_cast<value_type>(*this) + value;
        (*vector.vector)[0][index] = new_value;

        return *this;
      }



      inline const VectorReference &
      VectorReference::operator-=(const value_type &value) const
      {
        value_type new_value       = static_cast<value_type>(*this) - value;
        (*vector.vector)[0][index] = new_value;

        return *this;
      }



      inline const VectorReference &
      VectorReference::operator*=(const value_type &value) const
      {
        value_type new_value       = static_cast<value_type>(*this) * value;
        (*vector.vector)[0][index] = new_value;

        return *this;
      }



      inline const VectorReference &
      VectorReference::operator/=(const value_type &value) const
      {
        value_type new_value       = static_cast<value_type>(*this) / value;
        (*vector.vector)[0][index] = new_value;

        return *this;
      }
    } // namespace internal

#  endif /* DOXYGEN */


    inline internal::VectorReference
    Vector::operator()(const size_type index)
    {
      return internal::VectorReference(*this, index);
    }

    inline internal::VectorReference
    Vector::operator[](const size_type index)
    {
      return operator()(index);
    }

    inline Vector::value_type
    Vector::operator[](const size_type index) const
    {
      return operator()(index);
    }


    inline bool
    Vector::has_ghost_elements() const
    {
      return false;
    }
  } // namespace EpetraWrappers
} // namespace LinearAlgebra


/**
 * Declare dealii::LinearAlgebra::EpetraWrappers::Vector as distributed vector.
 */
template <>
struct is_serial_vector<LinearAlgebra::EpetraWrappers::Vector> : std::false_type
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
