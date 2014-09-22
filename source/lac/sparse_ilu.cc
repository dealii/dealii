// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#include <deal.II/lac/sparse_ilu.templates.h>

DEAL_II_NAMESPACE_OPEN


// explicit instantiations
template class SparseILU<double>;
template void SparseILU<double>::initialize<double> (const SparseMatrix<double> &,
                                                     const AdditionalData &data);
template void SparseILU<double>::decompose<double> (const SparseMatrix<double> &,
                                                    const double);
template void SparseILU<double>::vmult <double> (Vector<double> &,
                                                 const Vector<double> &) const;
template void SparseILU<double>::Tvmult <double> (Vector<double> &,
                                                  const Vector<double> &) const;
template void SparseILU<double>::initialize<float> (const SparseMatrix<float> &,
                                                    const AdditionalData &data);
template void SparseILU<double>::decompose<float> (const SparseMatrix<float> &,
                                                   const double);
template void SparseILU<double>::vmult<float> (Vector<float> &,
                                               const Vector<float> &) const;
template void SparseILU<double>::Tvmult<float> (Vector<float> &,
                                                const Vector<float> &) const;


template class SparseILU<float>;
template void SparseILU<float>::initialize<double> (const SparseMatrix<double> &,
                                                    const AdditionalData &data);
template void SparseILU<float>::decompose<double> (const SparseMatrix<double> &,
                                                   const double);
template void SparseILU<float>::vmult<double> (Vector<double> &,
                                               const Vector<double> &) const;
template void SparseILU<float>::Tvmult<double> (Vector<double> &,
                                                const Vector<double> &) const;
template void SparseILU<float>::initialize<float> (const SparseMatrix<float> &,
                                                   const AdditionalData &data);
template void SparseILU<float>::decompose<float> (const SparseMatrix<float> &,
                                                  const double);
template void SparseILU<float>::vmult<float> (Vector<float> &,
                                              const Vector<float> &) const;
template void SparseILU<float>::Tvmult<float> (Vector<float> &,
                                               const Vector<float> &) const;

DEAL_II_NAMESPACE_CLOSE
