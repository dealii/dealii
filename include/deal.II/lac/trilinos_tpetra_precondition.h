// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_precondition_h
#define dealii_trilinos_tpetra_precondition_h


#include <deal.II/base/config.h>

#include "deal.II/base/exceptions.h"
#include "deal.II/base/memory_space.h"

#include <Epetra_CombineMode.h>
#include <Teuchos_BLAS_types.hpp>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_TRILINOS_WITH_IFPACK2

#    include "deal.II/base/types.h"
#    include <deal.II/base/subscriptor.h>

#    include <deal.II/lac/la_parallel_vector.h>
#    include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#    include <Teuchos_ConfigDefs.hpp>
#    include <Teuchos_ParameterList.hpp>
#    include <Teuchos_RCPDecl.hpp>
#    include <Tpetra_Operator.hpp>


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    // forward declarations
#    ifndef DOXYGEN
    template <typename Number, typename MemorySpace>
    class SparseMatrix;
#    endif

    /**
     * The base class for all preconditioners based on Tpetra sparse matrices.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionBase : public Subscriptor
    {
    public:
      /**
       * Declare the type for container size.
       */
      using size_type = dealii::types::signed_global_dof_index;

      /**
       * Typedef for the NodeType
       */
#    if DEAL_II_TRILINOS_VERSION_GTE(14, 2, 0)
      using NodeType = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<
        typename MemorySpace::kokkos_space::execution_space,
        typename MemorySpace::kokkos_space>;
#    else
      using NodeType = Kokkos::Compat::KokkosDeviceWrapperNode<
        typename MemorySpace::kokkos_space::execution_space,
        typename MemorySpace::kokkos_space>;
#    endif

      /**
       * Standardized data struct to pipe additional flags to the
       * preconditioner.
       */
      struct AdditionalData
      {};

      /**
       * Constructor. Does not do anything. The <tt>initialize</tt> function of
       * the derived classes will have to create the preconditioner from a given
       * sparse matrix.
       */
      PreconditionBase() = default;

      /**
       * Destroys the preconditioner, leaving an object like just after having
       * called the constructor.
       */
      void
      clear();

      /**
       * Apply the preconditioner.
       */
      virtual void
      vmult(Vector<Number, MemorySpace>       &dst,
            const Vector<Number, MemorySpace> &src) const;

      /**
       * Apply the transpose preconditioner.
       */
      virtual void
      Tvmult(Vector<Number, MemorySpace>       &dst,
             const Vector<Number, MemorySpace> &src) const;

      /**
       * @name Access to underlying Trilinos data
       */
      /** @{ */
      /**
       *
       * Calling this function from an uninitialized object will cause an
       * exception.
       */
      Teuchos::RCP<Tpetra::Operator<Number, int, size_type, NodeType>>
      trilinos_operator() const;

      /** @} */

      /**
       * @name Partitioners
       */
      /** @{ */

      /**
       * Return the partitioning of the domain space of this matrix, i.e., the
       * partitioning of the vectors this matrix has to be multiplied with.
       */
      IndexSet
      locally_owned_domain_indices() const;

      /**
       * Return the partitioning of the range space of this matrix, i.e., the
       * partitioning of the vectors that are result from matrix-vector
       * products.
       */
      IndexSet
      locally_owned_range_indices() const;

      /** @} */

      /**
       * @addtogroup Exceptions
       */
      /** @{ */
      /**
       * Exception.
       */
      DeclException1(
        ExcNonMatchingMaps,
        std::string,
        << "The sparse matrix the preconditioner is based on "
        << "uses a map that is not compatible to the one in vector " << arg1
        << ". Check preconditioner and matrix setup.");

      /**
       * Exception.
       */
      DeclExceptionMsg(
        ExcTransposeNotSupported,
        "The chosen preconditioner does not support transposing the matrix.");
      /** @} */

    protected:
      /**
       * This is a pointer to the preconditioner object that is used when
       * applying the preconditioner.
       */
      Teuchos::RCP<Tpetra::Operator<Number, int, size_type, NodeType>>
        preconditioner;

      /**
       * An optional Teuchos::ParameterList for fine tuning the solver.
       * Please refer to the Amesos2 manual to see which parameters
       * to set for each individual solver.
       */
      Teuchos::ParameterList parameter_list;
    };

    /**
     * A wrapper class for a (pointwise) Jacobi preconditioner for Trilinos
     * Tpetra matrices. This preconditioner works both in serial and in parallel
     * depending on the matrix it is based on.
     *
     * The AdditionalData data structure allows to set preconditioner options.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionIfpackBase : public PreconditionBase<Number, MemorySpace>
    {
    public:
      PreconditionIfpackBase(const std::string &preconditioner_type);

      /**
       * Initializes the direct solver for the matrix <tt>A</tt> and creates a
       * factorization for it with the package chosen from the additional
       * data structure. Note that there is no need for a preconditioner
       * here and solve() is not called.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A);

    protected:
      /**
       *  The set preconditioner type to be handed to the factory of Ifpack2
       */
      std::string preconditioner_type;
    };


    /**
     * The class for the classical Jacobi preconditioner within Ifpack2.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionJacobi
      : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      struct AdditionalData
      {
        AdditionalData(const double omega        = 1,
                       const int    n_sweeps     = 1,
                       const bool   fix_diagonal = false,
                       const double min_diagonal = 0);

        double omega;
        int    n_sweeps;
        bool   fix_diagonal;
        double min_diagonal;
      };

      PreconditionJacobi();

      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };

    /**
     * The class for the l1 variant of the Jacobi preconditioner within Ifpack2,
     * first described in TODO: citation Saad et. al.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionL1Jacobi
      : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      struct AdditionalData
      {
        AdditionalData(const double omega        = 1,
                       const int    n_sweeps     = 1,
                       const bool   fix_diagonal = false,
                       const double min_diagonal = 0);

        double omega;
        int    n_sweeps;
        bool   fix_diagonal;
        double min_diagonal;
      };

      PreconditionL1Jacobi();

      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };

  } // namespace TpetraWrappers
} // namespace LinearAlgebra


/** @} */


DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_TRILINOS_WITH_IFPACK2
#endif   // DEAL_II_TRILINOS_WITH_TPETRA

#endif
