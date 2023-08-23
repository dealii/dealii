// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2023 by the deal.II authors
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

#ifndef dealii_trilinos_epetra_vector_h
#define dealii_trilinos_epetra_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/mpi_stub.h>
#  include <deal.II/base/subscriptor.h>

#  include <deal.II/lac/read_vector.h>
#  include <deal.II/lac/trilinos_epetra_communication_pattern.h>
#  include <deal.II/lac/vector_operation.h>
#  include <deal.II/lac/vector_type_traits.h>

#  include <Epetra_FEVector.h>

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
     * This class implements a wrapper to the Trilinos distributed vector
     * class Epetra_FEVector. This class requires Trilinos to be compiled
     * with MPI support.
     *
     * @ingroup TrilinosWrappers
     * @ingroup Vectors
     */
    class Vector : public ReadVector<double>, public Subscriptor
    {
    public:
      using value_type = double;
      using size_type  = types::global_dof_index;
      using real_type  = double;

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
       * filled with zero (false) or left untouched (true).
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
        ArrayView<double> &elements) const override;

      /**
       * Copy function. This function takes a Vector and copies all the
       * elements. The Vector will have the same parallel distribution as @p
       * V.
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
    };


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

#endif

#endif
