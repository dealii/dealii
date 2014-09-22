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

#include <deal.II/lac/sparse_vanka.templates.h>

DEAL_II_NAMESPACE_OPEN


// explicit instantiations
template class SparseVanka<float>;
template class SparseVanka<double>;

template void SparseVanka<double>::vmult<float> (Vector<float>       &dst,
                                                 const Vector<float> &src) const;
template void SparseVanka<double>::vmult<double> (Vector<double>       &dst,
                                                  const Vector<double> &src) const;


template class SparseBlockVanka<float>;
template class SparseBlockVanka<double>;

template void SparseBlockVanka<double>::vmult<float> (Vector<float>       &dst,
                                                      const Vector<float> &src) const;
template void SparseBlockVanka<double>::vmult<double> (Vector<double>       &dst,
                                                       const Vector<double> &src) const;

DEAL_II_NAMESPACE_CLOSE
