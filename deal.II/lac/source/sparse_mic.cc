//----------------------------  sparse_mic.cc  ---------------------------
//    Copyright (C) 1998, 1999, 2000, 2001, 2002
//    by the deal.II authors and Stephen "Cheffo" Kolaroff
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_mic.cc  ---------------------------

#include <lac/sparse_mic.templates.h>


// explicit instantiations for double and float matrices
template class SparseMIC<double>;
template void SparseMIC<double>::decompose<double> (const SparseMatrix<double> &,
                                                    const double);
template void SparseMIC<double>::vmult<double> (Vector<double> &,
                                                const Vector<double> &) const;
template void SparseMIC<double>::decompose<float> (const SparseMatrix<float> &,
                                                   const double);
template void SparseMIC<double>::vmult<float> (Vector<float> &,
                                               const Vector<float> &) const;

template class SparseMIC<float>;
template void SparseMIC<float>::decompose<double> (const SparseMatrix<double> &,
                                                   const double);
template void SparseMIC<float>::vmult<double> (Vector<double> &,
                                               const Vector<double> &) const;
template void SparseMIC<float>::decompose<float> (const SparseMatrix<float> &,
                                                  const double);
template void SparseMIC<float>::vmult<float> (Vector<float> &,
                                              const Vector<float> &) const;

