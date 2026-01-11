// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_to_trilinos_wrappers_h
#define dealii_trilinos_tpetra_to_trilinos_wrappers_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_TRILINOS

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
namespace LinearAlgebra::TpetraWrappers
{
  template <typename Number, typename MemorySpace>
  class BlockSparseMatrix;

  template <typename Number, typename MemorySpace>
  class SparseMatrix;

  template <typename MemorySpace>
  class SparsityPattern;

  template <typename Number, typename MemorySpace>
  class BlockVector;

  template <typename Number, typename MemorySpace>
  class Vector;

  template <typename Number, typename MemorySpace>
  class PreconditionBase;

  template <typename Number, typename MemorySpace>
  class PreconditionIdentity;

  template <typename Number, typename MemorySpace>
  class PreconditionIfpackBase;

  template <typename Number, typename MemorySpace>
  class PreconditionIfpack;

  template <typename Number, typename MemorySpace>
  class PreconditionJacobi;

  template <typename Number, typename MemorySpace>
  class PreconditionL1Jacobi;

  template <typename Number, typename MemorySpace>
  class PreconditionL1GaussSeidel;

  template <typename Number, typename MemorySpace>
  class PreconditionSOR;

  template <typename Number, typename MemorySpace>
  class PreconditionSSOR;

  template <typename Number, typename MemorySpace>
  class PreconditionChebyshev;

  template <typename Number, typename MemorySpace>
  class PreconditionILU;

  template <typename Number, typename MemorySpace>
  class PreconditionILUT;

  template <typename Number, typename MemorySpace>
  class PreconditionBlockJacobi;

  template <typename Number, typename MemorySpace>
  class PreconditionBlockSOR;

  template <typename Number, typename MemorySpace>
  class PreconditionBlockSSOR;
} // namespace LinearAlgebra::TpetraWrappers
#  endif

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
namespace MemorySpace
{
  struct Host;
}
#  endif

namespace TrilinosWrappers
{
  namespace MPI
  {
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
    using Vector = ::dealii::LinearAlgebra::TpetraWrappers::
      Vector<double, ::dealii::MemorySpace::Host>;
    using BlockVector = ::dealii::LinearAlgebra::TpetraWrappers::
      BlockVector<double, ::dealii::MemorySpace::Host>;
#  else
    class Vector;
    class BlockVector;
#  endif
  } // namespace MPI

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
  using BlockSparseMatrix = ::dealii::LinearAlgebra::TpetraWrappers::
    BlockSparseMatrix<double, ::dealii::MemorySpace::Host>;
  using SparseMatrix = ::dealii::LinearAlgebra::TpetraWrappers::
    SparseMatrix<double, ::dealii::MemorySpace::Host>;
  using SparsityPattern =
    ::dealii::LinearAlgebra::TpetraWrappers::SparsityPattern<
      ::dealii::MemorySpace::Host>;
#  else
  class BlockSparseMatrix;
  class SparseMatrix;
  class SparsityPattern;
#  endif

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
  using PreconditionBase = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionBase<double, MemorySpace::Host>;
  using PreconditionIdentity = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionIdentity<double, MemorySpace::Host>;
  using PreconditionIfpackBase = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionIfpackBase<double, MemorySpace::Host>;
  using PreconditionIfpack = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionIfpack<double, MemorySpace::Host>;
  using PreconditionJacobi = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionJacobi<double, MemorySpace::Host>;
  using PreconditionL1Jacobi = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionL1Jacobi<double, MemorySpace::Host>;
  using PreconditionL1GaussSeidel = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionL1GaussSeidel<double, MemorySpace::Host>;
  using PreconditionSOR =
    ::dealii::LinearAlgebra::TpetraWrappers::PreconditionSOR<double,
                                                             MemorySpace::Host>;
  using PreconditionSSOR = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionSSOR<double, MemorySpace::Host>;
  using PreconditionChebyshev = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionChebyshev<double, MemorySpace::Host>;
  using PreconditionILU =
    ::dealii::LinearAlgebra::TpetraWrappers::PreconditionILU<double,
                                                             MemorySpace::Host>;
  using PreconditionILUT = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionILUT<double, MemorySpace::Host>;
  using PreconditionBlockJacobi = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionBlockJacobi<double, MemorySpace::Host>;
  using PreconditionBlockSOR = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionBlockSOR<double, MemorySpace::Host>;
  using PreconditionBlockSSOR = ::dealii::LinearAlgebra::TpetraWrappers::
    PreconditionBlockSSOR<double, MemorySpace::Host>;
#  else
  class PreconditionBase;
  class PreconditionIdentity;
  class PreconditionIfpackBase;
  class PreconditionIfpack;
  class PreconditionJacobi;
  class PreconditionL1Jacobi;
  class PreconditionL1GaussSeidel;
  class PreconditionSOR;
  class PreconditionSSOR;
  class PreconditionChebyshev;
  class PreconditionILU;
  class PreconditionILUT;
  class PreconditionBlockJacobi;
  class PreconditionBlockSOR;
  class PreconditionBlockSSOR;
#  endif
} // namespace TrilinosWrappers

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
