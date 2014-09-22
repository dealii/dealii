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


#include <deal.II/lac/sparse_mic.templates.h>

DEAL_II_NAMESPACE_OPEN


// explicit instantiations for double and float matrices
template class SparseMIC<double>;
template void SparseMIC<double>::initialize<double> (const SparseMatrix<double> &,
                                                     const AdditionalData &data);
template void SparseMIC<double>::decompose<double> (const SparseMatrix<double> &,
                                                    const double);
template void SparseMIC<double>::vmult<double> (Vector<double> &,
                                                const Vector<double> &) const;
template void SparseMIC<double>::initialize<float> (const SparseMatrix<float> &,
                                                    const AdditionalData &data);
template void SparseMIC<double>::decompose<float> (const SparseMatrix<float> &,
                                                   const double);
template void SparseMIC<double>::vmult<float> (Vector<float> &,
                                               const Vector<float> &) const;

template class SparseMIC<float>;
template void SparseMIC<float>::initialize<double> (const SparseMatrix<double> &,
                                                    const AdditionalData &data);
template void SparseMIC<float>::decompose<double> (const SparseMatrix<double> &,
                                                   const double);
template void SparseMIC<float>::vmult<double> (Vector<double> &,
                                               const Vector<double> &) const;
template void SparseMIC<float>::initialize<float> (const SparseMatrix<float> &,
                                                   const AdditionalData &data);
template void SparseMIC<float>::decompose<float> (const SparseMatrix<float> &,
                                                  const double);
template void SparseMIC<float>::vmult<float> (Vector<float> &,
                                              const Vector<float> &) const;


DEAL_II_NAMESPACE_CLOSE
