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

#ifndef dealii_trilinos_tpetra_types_h
#define dealii_trilinos_tpetra_types_h

#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/base/memory_space.h>

#  include <Kokkos_DualView.hpp>
#  include <Tpetra_Details_DefaultTypes.hpp>


// Forward declarations
#  ifndef DOXYGEN
#    ifdef DEAL_II_TRILINOS_WITH_IFPACK2
#      include <Ifpack2_Preconditioner.hpp>
#    endif
#    include <Tpetra_CrsGraph_fwd.hpp>
#    include <Tpetra_CrsMatrix_fwd.hpp>
#    include <Tpetra_Export_fwd.hpp>
#    include <Tpetra_Import_fwd.hpp>
#    include <Tpetra_Map_fwd.hpp>
#    include <Tpetra_MultiVector_fwd.hpp>
#    include <Tpetra_Operator_fwd.hpp>
#    include <Tpetra_RowMatrix_fwd.hpp>
#    include <Tpetra_Vector_fwd.hpp>
#  endif

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {

    namespace TpetraTypes
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
      //
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
       * Communication between processors.
       */
      template <typename MemorySpace>
      using ExportType = Tpetra::Export<LO, GO, NodeType<MemorySpace>>;

      /**
       * Communication between processors.
       */
      template <typename MemorySpace>
      using ImportType = Tpetra::Import<LO, GO, NodeType<MemorySpace>>;

      /**
       * Tpetra equivalent of IndexSet.
       */
      template <typename MemorySpace>
      using MapType = Tpetra::Map<LO, GO, NodeType<MemorySpace>>;

      /**
       * Tpetra sparsity pattern type.
       */
      template <typename MemorySpace>
      using GraphType = Tpetra::CrsGraph<LO, GO, NodeType<MemorySpace>>;

      /**
       * Tpetra Vector type.
       */
      template <typename Number, typename MemorySpace>
      using VectorType = Tpetra::Vector<Number, LO, GO, NodeType<MemorySpace>>;

      /**
       * Tpetra type for a row column vectors.
       */
      template <typename Number, typename MemorySpace>
      using MultiVectorType =
        Tpetra::MultiVector<Number, LO, GO, NodeType<MemorySpace>>;


      /**
       * General Tpetra class for a linear operator,
       * e.g. a Matrix or Preconditioner.
       *
       */
      template <typename Number, typename MemorySpace>
      using LinearOperator =
        Tpetra::Operator<Number, LO, GO, NodeType<MemorySpace>>;


      /**
       * Tpetra type for a parallel distributed sparse matrix in Crs or CSR
       * format.
       */
      template <typename Number, typename MemorySpace>
      using MatrixType =
        Tpetra::CrsMatrix<Number, LO, GO, NodeType<MemorySpace>>;


      /**
       * Tpetra type for a parallel distributed row matrix.
       */
      template <typename Number, typename MemorySpace>
      using RowMatrixType =
        Tpetra::RowMatrix<Number, LO, GO, NodeType<MemorySpace>>;

      /**
       * @brief Typedef for the Kokkos::DualView type.
       * This is needed for shallow copies of deal.II LA structures
       * to Trilinos LA structures.
       */
      template <typename Number, typename MemorySpace>
      using DualViewType =
        typename VectorType<Number, MemorySpace>::dual_view_type;

      /**
       * @brief Typedef for the Kokkos::View type.
       * This is needed for shallow copies of deal.II LA structures
       * to Trilinos LA structures.
       *
       */
      template <typename Number>
      using HostViewType =
        typename DualViewType<Number, dealii::MemorySpace::Host>::t_host;

#  ifdef DEAL_II_TRILINOS_WITH_IFPACK2
      /**
       * Type for a Trilinos preconditioner from the Ifpack2 package.
       */
      template <typename Number, typename MemorySpace>
      using Ifpack2PreconType =
        Ifpack2::Preconditioner<Number, LO, GO, NodeType<MemorySpace>>;
#  endif

    } // namespace TpetraTypes
  }   // namespace TpetraWrappers
} // namespace LinearAlgebra
DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
