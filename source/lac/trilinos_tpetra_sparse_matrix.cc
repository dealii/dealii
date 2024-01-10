// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/lac/trilinos_tpetra_sparse_matrix.templates.h>

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN
// explicit instantiations
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    template class SparseMatrix<double>;

    template void
    SparseMatrix<double>::reinit(
      const IndexSet                       &row_parallel_partitioning,
      const IndexSet                       &col_parallel_partitioning,
      const dealii::DynamicSparsityPattern &sparsity_pattern,
      const MPI_Comm                        communicator,
      const bool                            exchange_data);

  } // namespace TpetraWrappers
} // namespace LinearAlgebra
#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA
