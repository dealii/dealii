// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/lac/sparse_decomposition.templates.h>

DEAL_II_NAMESPACE_OPEN


template class SparseLUDecomposition<double>;
template void SparseLUDecomposition<double>::initialize<double> (const SparseMatrix<double> &,
    const AdditionalData data);
template void SparseLUDecomposition<double>::initialize<float> (const SparseMatrix<float> &,
    const AdditionalData data);
template void SparseLUDecomposition<double>::decompose<double> (const SparseMatrix<double> &,
    const double);
template void SparseLUDecomposition<double>::decompose<float> (const SparseMatrix<float> &,
    const double);

template void SparseLUDecomposition<double>::copy_from<double> (const SparseMatrix<double> &);
template void SparseLUDecomposition<double>::copy_from<float> (const SparseMatrix<float> &);


template class SparseLUDecomposition<float>;
template void SparseLUDecomposition<float>::initialize<double> (const SparseMatrix<double> &,
    const AdditionalData data);
template void SparseLUDecomposition<float>::initialize<float> (const SparseMatrix<float> &,
    const AdditionalData data);
template void SparseLUDecomposition<float>::decompose<double> (const SparseMatrix<double> &,
    const double);
template void SparseLUDecomposition<float>::decompose<float> (const SparseMatrix<float> &,
    const double);

template void SparseLUDecomposition<float>::copy_from<double> (const SparseMatrix<double> &);
template void SparseLUDecomposition<float>::copy_from<float> (const SparseMatrix<float> &);

DEAL_II_NAMESPACE_CLOSE
