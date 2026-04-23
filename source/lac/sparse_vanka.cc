// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 1999 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/sparse_vanka.h>
#include <deal.II/lac/sparse_vanka.templates.h>

#include <taskflow/core/async.hpp>
#include <taskflow/core/graph.hpp>
#include <taskflow/utility/traits.hpp>

#include <ostream>

DEAL_II_NAMESPACE_OPEN


// explicit instantiations
template class SparseVanka<float>;
template class SparseVanka<double>;

template void
SparseVanka<double>::vmult<float>(Vector<float>       &dst,
                                  const Vector<float> &src) const;
template void
SparseVanka<double>::vmult<double>(Vector<double>       &dst,
                                   const Vector<double> &src) const;
template void
SparseVanka<double>::Tvmult<float>(Vector<float>       &dst,
                                   const Vector<float> &src) const;
template void
SparseVanka<double>::Tvmult<double>(Vector<double>       &dst,
                                    const Vector<double> &src) const;


template class SparseBlockVanka<float>;
template class SparseBlockVanka<double>;

template void
SparseBlockVanka<double>::vmult<float>(Vector<float>       &dst,
                                       const Vector<float> &src) const;
template void
SparseBlockVanka<double>::vmult<double>(Vector<double>       &dst,
                                        const Vector<double> &src) const;
template void
SparseBlockVanka<double>::Tvmult<float>(Vector<float>       &dst,
                                        const Vector<float> &src) const;
template void
SparseBlockVanka<double>::Tvmult<double>(Vector<double>       &dst,
                                         const Vector<double> &src) const;

DEAL_II_NAMESPACE_CLOSE
