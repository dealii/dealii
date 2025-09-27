// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/sparse_ilu.templates.h>

DEAL_II_NAMESPACE_OPEN


// explicit instantiations
template class SparseILU<double>;
template void
SparseILU<double>::initialize<double>(const SparseMatrix<double> &,
                                      const AdditionalData &data);
template void
SparseILU<double>::vmult<double>(Vector<double> &,
                                 const Vector<double> &) const;
template void
SparseILU<double>::Tvmult<double>(Vector<double> &,
                                  const Vector<double> &) const;
template void
SparseILU<double>::initialize<float>(const SparseMatrix<float> &,
                                     const AdditionalData &data);
template void
SparseILU<double>::vmult<float>(Vector<float> &, const Vector<float> &) const;
template void
SparseILU<double>::Tvmult<float>(Vector<float> &, const Vector<float> &) const;


template class SparseILU<float>;
template void
SparseILU<float>::initialize<double>(const SparseMatrix<double> &,
                                     const AdditionalData &data);
template void
SparseILU<float>::vmult<double>(Vector<double> &, const Vector<double> &) const;
template void
SparseILU<float>::Tvmult<double>(Vector<double> &,
                                 const Vector<double> &) const;
template void
SparseILU<float>::initialize<float>(const SparseMatrix<float> &,
                                    const AdditionalData &data);
template void
SparseILU<float>::vmult<float>(Vector<float> &, const Vector<float> &) const;
template void
SparseILU<float>::Tvmult<float>(Vector<float> &, const Vector<float> &) const;

DEAL_II_NAMESPACE_CLOSE
