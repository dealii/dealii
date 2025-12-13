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

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
namespace LinearAlgebra::TpetraWrappers
{
  template <typename Number, typename MemorySpace>
  class SparseMatrix;

  template <typename MemorySpace>
  class SparsityPattern;

  template <typename Number, typename MemorySpace>
  class Vector;

  template <typename Number, typename MemorySpace>
  class BlockVector;
} // namespace LinearAlgebra::TpetraWrappers
#endif

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
namespace MemorySpace
{
  struct Host;
}
#endif

namespace TrilinosWrappers
{
  namespace MPI
  {
#ifdef DEAL_II_TRILINOS_WITH_TPETRA
    using SparseMatrix = ::dealii::LinearAlgebra::TpetraWrappers::
      SparseMatrix<double, ::dealii::MemorySpace::Host>;

    using Vector = ::dealii::LinearAlgebra::TpetraWrappers::
      Vector<double, ::dealii::MemorySpace::Host>;
    using BlockVector = ::dealii::LinearAlgebra::TpetraWrappers::
      BlockVector<double, ::dealii::MemorySpace::Host>;
#else
    class SparseMatrix;
    class Vector;
    class BlockVector;
#endif
  } // namespace MPI
} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
