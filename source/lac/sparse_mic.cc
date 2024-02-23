// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/sparse_mic.templates.h>

DEAL_II_NAMESPACE_OPEN


// explicit instantiations for double and float matrices
template class SparseMIC<double>;
template void
SparseMIC<double>::initialize<double>(const SparseMatrix<double> &,
                                      const AdditionalData &data);
template void
SparseMIC<double>::vmult<double>(Vector<double> &,
                                 const Vector<double> &) const;
template void
SparseMIC<double>::Tvmult<double>(Vector<double> &,
                                  const Vector<double> &) const;
template void
SparseMIC<double>::initialize<float>(const SparseMatrix<float> &,
                                     const AdditionalData &data);
template void
SparseMIC<double>::vmult<float>(Vector<float> &, const Vector<float> &) const;
template void
SparseMIC<double>::Tvmult<float>(Vector<float> &, const Vector<float> &) const;

template class SparseMIC<float>;
template void
SparseMIC<float>::initialize<double>(const SparseMatrix<double> &,
                                     const AdditionalData &data);
template void
SparseMIC<float>::vmult<double>(Vector<double> &, const Vector<double> &) const;
template void
SparseMIC<float>::Tvmult<double>(Vector<double> &,
                                 const Vector<double> &) const;
template void
SparseMIC<float>::initialize<float>(const SparseMatrix<float> &,
                                    const AdditionalData &data);
template void
SparseMIC<float>::vmult<float>(Vector<float> &, const Vector<float> &) const;
template void
SparseMIC<float>::Tvmult<float>(Vector<float> &, const Vector<float> &) const;



DEAL_II_NAMESPACE_CLOSE
