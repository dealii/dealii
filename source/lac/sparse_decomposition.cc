// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/lac/sparse_decomposition.templates.h>

DEAL_II_NAMESPACE_OPEN


template class SparseLUDecomposition<double>;
template void
SparseLUDecomposition<double>::initialize<double>(const SparseMatrix<double> &,
                                                  const AdditionalData data);
template void
SparseLUDecomposition<double>::initialize<float>(const SparseMatrix<float> &,
                                                 const AdditionalData data);

template void
SparseLUDecomposition<double>::copy_from<double>(const SparseMatrix<double> &);
template void
SparseLUDecomposition<double>::copy_from<float>(const SparseMatrix<float> &);


template class SparseLUDecomposition<float>;
template void
SparseLUDecomposition<float>::initialize<double>(const SparseMatrix<double> &,
                                                 const AdditionalData data);
template void
SparseLUDecomposition<float>::initialize<float>(const SparseMatrix<float> &,
                                                const AdditionalData data);

template void
SparseLUDecomposition<float>::copy_from<double>(const SparseMatrix<double> &);
template void
SparseLUDecomposition<float>::copy_from<float>(const SparseMatrix<float> &);

DEAL_II_NAMESPACE_CLOSE
