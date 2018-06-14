// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
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

#include <deal.II/lac/matrix_lib.templates.h>
#include <deal.II/lac/sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

MeanValueFilter::MeanValueFilter(size_type component)
  : component(component)
{}

template void
MeanValueFilter::filter(Vector<float> &) const;
template void
MeanValueFilter::filter(Vector<double> &) const;
template void
MeanValueFilter::filter(BlockVector<float> &) const;
template void
MeanValueFilter::filter(BlockVector<double> &) const;
template void
MeanValueFilter::vmult(Vector<float> &, const Vector<float> &) const;
template void
MeanValueFilter::vmult(Vector<double> &, const Vector<double> &) const;
template void
MeanValueFilter::vmult(BlockVector<float> &, const BlockVector<float> &) const;
template void
MeanValueFilter::vmult(BlockVector<double> &,
                       const BlockVector<double> &) const;

template void
MeanValueFilter::vmult_add(Vector<float> &, const Vector<float> &) const;
template void
MeanValueFilter::vmult_add(Vector<double> &, const Vector<double> &) const;
template void
MeanValueFilter::vmult_add(BlockVector<float> &,
                           const BlockVector<float> &) const;
template void
MeanValueFilter::vmult_add(BlockVector<double> &,
                           const BlockVector<double> &) const;

DEAL_II_NAMESPACE_CLOSE
