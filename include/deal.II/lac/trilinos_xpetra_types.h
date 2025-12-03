// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_xpetra_types_h
#define dealii_trilinos_xpetra_types_h

#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#ifdef DEAL_II_TRILINOS_WITH_XPETRA
// Forward declarations
#  ifndef DOXYGEN
#    include <Xpetra_CrsGraphFactory.hpp>
#    include <Xpetra_CrsMatrixWrap.hpp>
#    include <Xpetra_DefaultPlatform.hpp>
#    include <Xpetra_Map.hpp>
#    include <Xpetra_MultiVectorFactory_decl.hpp>
#    include <Xpetra_Parameters.hpp>

#    ifdef DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH
//     //FROSch
#      include <FROSch_OneLevelPreconditioner_def.hpp>
#      include <FROSch_SchwarzPreconditioners_fwd.hpp>
#      include <FROSch_Tools_def.hpp>
#      include <FROSch_TwoLevelPreconditioner_def.hpp>
#      include <ShyLU_DDFROSch_config.h>
#    endif // DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH
#  endif   // DOXYGEN



DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {

    namespace XpetraTypes
    {
      /**
       * local ordinate (processor local indices).
       */
      using LO = int;

      /**
       * global ordinate (global indices).
       */
      using GO = types::signed_global_dof_index;

#  if DEAL_II_TRILINOS_VERSION_GTE(14, 2, 0)
      /**
       * Where and how calculations should be executed,
       * i.e. Host (Serial,OpenMP) or Device (GPU)
       */
      template <typename MemorySpace>
      using NodeType = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<
        typename MemorySpace::kokkos_space::execution_space,
        typename MemorySpace::kokkos_space>;
#  else
      template <typename MemorySpace>
      using NodeType = Kokkos::Compat::KokkosDeviceWrapperNode<
        typename MemorySpace::kokkos_space::execution_space,
        typename MemorySpace::kokkos_space>;
#  endif

      /**
       * Xpetra equivalent of IndexSet.
       */
      template <typename MemorySpace>
      using MapType = Xpetra::Map<LO, GO, NodeType<MemorySpace>>;

      /**
       * Xpetra sparsity pattern type.
       */
      template <typename MemorySpace>
      using GraphType = Xpetra::CrsGraph<LO, GO, NodeType<MemorySpace>>;

      /**
       * Xpetra MultiVector type.
       */
      template <typename Number, typename MemorySpace>
      using MultiVectorType =
        Xpetra::MultiVector<Number, LO, GO, NodeType<MemorySpace>>;

      /**
       * General Xpetra class for a linear operator,
       * e.g. a matrix or preconditioner.
       */
      template <typename Number, typename MemorySpace>
      using LinearOperator =
        Xpetra::Operator<Number, LO, GO, NodeType<MemorySpace>>;

      /**
       * Xpetra type for a parallel distributed sparse matrix.
       */
      template <typename Number, typename MemorySpace>
      using MatrixType = Xpetra::Matrix<Number, LO, GO, NodeType<MemorySpace>>;

      /**
       * Xpetra type for a parallel distributed sparse matrix in CSR format.
       */
      template <typename Number, typename MemorySpace>
      using CrsMatrixType =
        Xpetra::CrsMatrix<Number, LO, GO, NodeType<MemorySpace>>;

      /**
       * Factory to create Xpetra::Maps
       */
      template <typename MemorySpace>
      using MapFactoryType = Xpetra::MapFactory<LO, GO, NodeType<MemorySpace>>;

      /**
       * Factory to create Xpetra::CrsGraphs.
       */
      template <typename MemorySpace>
      using GraphFactoryType =
        Xpetra::CrsGraphFactory<LO, GO, NodeType<MemorySpace>>;

      /**
       * Factory to create Xpetra::Multivectors.
       */
      template <typename Number, typename MemorySpace>
      using MultiVectorFactoryType =
        Xpetra::MultiVectorFactory<Number, LO, GO, NodeType<MemorySpace>>;

      /**
       * Xpetra type to transfer a Tpetra::Map into a Xpetra::Map
       */
      template <typename MemorySpace>
      using TpetraMapType = Xpetra::TpetraMap<LO, GO, NodeType<MemorySpace>>;

      /**
       * Xpetra type to transfer a Tpetra::CrsGraph into a Xpetra::CrsGraph
       */
      template <typename MemorySpace>
      using TpetraGraphType =
        Xpetra::TpetraCrsGraph<LO, GO, NodeType<MemorySpace>>;

      /**
       * Xpetra type to transfer a Tpetra::MultiVector into a
       * Xpetra::MultiVector
       */
      template <typename Number, typename MemorySpace>
      using XTpetraMultiVectorType =
        Xpetra::TpetraMultiVector<Number, LO, GO, NodeType<MemorySpace>>;

      /**
       * Xpetra type to transfer a Tpetra::CrsMatrix into an Xpetra::CrsMatrix
       */
      template <typename Number, typename MemorySpace>
      using TpetraCrsMatrixType =
        Xpetra::TpetraCrsMatrix<Number, LO, GO, NodeType<MemorySpace>>;

      /**
       * Xpetra helper type to transfer a Tpetra::CrsMatrix into an
       * Xpetra::CrsMatrix
       */
      template <typename Number, typename MemorySpace>
      using CrsMatrixWrapType =
        Xpetra::CrsMatrixWrap<Number, LO, GO, NodeType<MemorySpace>>;

#  ifdef DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH
      /**
       * Adapter for Xpetra::Operator to Tpetra::Operator.
       */
      template <typename Number, typename MemorySpace>
      using XpetraToTpetraLinearOperator =
        FROSch::TpetraPreconditioner<Number, LO, GO, NodeType<MemorySpace>>;

      /**
       * Type for a Trilinos preconditioner from the FROSch package.
       */
      template <typename Number, typename MemorySpace>
      using FROSchOneLevelType =
        FROSch::OneLevelPreconditioner<Number, LO, GO, NodeType<MemorySpace>>;

      template <typename Number, typename MemorySpace>
      using FROSchTwoLevelType =
        FROSch::TwoLevelPreconditioner<Number, LO, GO, NodeType<MemorySpace>>;
#  endif // DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH

    } // namespace XpetraTypes
  }   // namespace TpetraWrappers
} // namespace LinearAlgebra
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_XPETRA

#endif // dealii_trilinos_xpetra_types_h
