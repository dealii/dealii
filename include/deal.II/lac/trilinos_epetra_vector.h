// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2019 by the deal.II authors
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

#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/subscriptor.h>

#  include <deal.II/lac/trilinos_epetra_communication_pattern.h>
#  include <deal.II/lac/vector_operation.h>
#  include <deal.II/lac/vector_space_vector.h>
#  include <deal.II/lac/vector_type_traits.h>

#  include <Epetra_FEVector.h>
#  include <mpi.h>

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
     * class Epetra_FEVector. This class is derived from the
     * LinearAlgebra::VectorSpaceVector class. Note however that Epetra only
     * works with Number = double. This class requires Trilinos to be compiled
     * with MPI support.
     *
     * @ingroup TrilinosWrappers
     * @ingroup Vectors
     * @author Bruno Turcksin, 2015
     */
    class Vector : public VectorSpaceVector<double>, public Subscriptor
    {
    public:
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
                      const MPI_Comm &communicator);

      /**
       * Reinit functionality. This function destroys the old vector content
       * and generates a new one based on the input partitioning. The flag
       * <tt>omit_zeroing_entries</tt> determines whether the vector should be
       * filled with zero (false) or left untouched (true).
       */
      void
      reinit(const IndexSet &parallel_partitioner,
             const MPI_Comm &communicator,
             const bool      omit_zeroing_entries = false);

      /**
       * Change the dimension to that of the vector @p V. The elements of @p V are not
       * copied.
       */
      virtual void
      reinit(const VectorSpaceVector<double> &V,
             const bool omit_zeroing_entries = false) override;

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
      virtual Vector &
      operator=(const double s) override;

      /**
       * Imports all the elements present in the vector's IndexSet from the
       * input vector @p V. VectorOperation::values @p operation is used to decide if
       * the elements in @p V should be added to the current vector or replace the
       * current elements. The last parameter can be used if the same
       * communication pattern is used multiple times. This can be used to
       * improve performance.
       */
      virtual void
      import(
        const ReadWriteVector<double> &                 V,
        VectorOperation::values                         operation,
        std::shared_ptr<const CommunicationPatternBase> communication_pattern =
          std::shared_ptr<const CommunicationPatternBase>()) override;

      /**
       * Multiply the entire vector by a fixed factor.
       */
      virtual Vector &
      operator*=(const double factor) override;

      /**
       * Divide the entire vector by a fixed factor.
       */
      virtual Vector &
      operator/=(const double factor) override;

      /**
       * Add the vector @p V to the present one.
       */
      virtual Vector &
      operator+=(const VectorSpaceVector<double> &V) override;

      /**
       * Subtract the vector @p V from the present one.
       */
      virtual Vector &
      operator-=(const VectorSpaceVector<double> &V) override;

      /**
       * Return the scalar product of two vectors. The vectors need to have the
       * same layout.
       */
      virtual double
      operator*(const VectorSpaceVector<double> &V) const override;

      /**
       * Add @p a to all components. Note that @p is a scalar not a vector.
       */
      virtual void
      add(const double a) override;

      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this +=
       * a*V</tt>. The vectors need to have the same layout.
       */
      virtual void
      add(const double a, const VectorSpaceVector<double> &V) override;

      /**
       * Multiple addition of multiple of a vector, i.e. <tt>*this> +=
       * a*V+b*W</tt>. The vectors need to have the same layout.
       */
      virtual void
      add(const double                     a,
          const VectorSpaceVector<double> &V,
          const double                     b,
          const VectorSpaceVector<double> &W) override;

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this
       * = s*(*this)+a*V</tt>.
       */
      virtual void
      sadd(const double                     s,
           const double                     a,
           const VectorSpaceVector<double> &V) override;

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument. This function is mostly meant to simulate multiplication
       * (and immediate re-assignment) by a diagonal scaling matrix. The
       * vectors need to have the same layout.
       */
      virtual void
      scale(const VectorSpaceVector<double> &scaling_factors) override;

      /**
       * Assignment <tt>*this = a*V</tt>.
       */
      virtual void
      equ(const double a, const VectorSpaceVector<double> &V) override;

      /**
       * Return whether the vector contains only elements with value zero.
       */
      virtual bool
      all_zero() const override;

      /**
       * Return the mean value of the element of this vector.
       */
      virtual double
      mean_value() const override;

      /**
       * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
       * absolute values of all entries among all processors).
       */
      virtual double
      l1_norm() const override;

      /**
       * Return the l<sub>2</sub> norm of the vector (i.e., the square root of
       * the sum of the square of all entries among all processors).
       */
      virtual double
      l2_norm() const override;

      /**
       * Return the maximum norm of the vector (i.e., the maximum absolute value
       * among all entries and among all processors).
       */
      virtual double
      linfty_norm() const override;

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
      virtual double
      add_and_dot(const double                     a,
                  const VectorSpaceVector<double> &V,
                  const VectorSpaceVector<double> &W) override;
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
       * Return the MPI communicator object in use with this object.
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
      virtual ::dealii::IndexSet
      locally_owned_elements() const override;

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
                                 const MPI_Comm &mpi_comm);

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
