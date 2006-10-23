//---------------------------------------------------------------------------
//    Copyright (C) 2002, 2003, 2005, 2006 by the deal.II authors
//    by the deal.II authors and Stephen "Cheffo" Kolaroff
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/sparse_mic.templates.h>

DEAL_II_NAMESPACE_OPEN


// explicit instantiations for double and float matrices
template class SparseMIC<double>;
template void SparseMIC<double>::initialize<double> (const SparseMatrix<double> &,
						     const AdditionalData data);
template void SparseMIC<double>::decompose<double> (const SparseMatrix<double> &,
                                                    const double);
template void SparseMIC<double>::vmult<double> (Vector<double> &,
                                                const Vector<double> &) const;
template void SparseMIC<double>::initialize<float> (const SparseMatrix<float> &,
						    const AdditionalData data);
template void SparseMIC<double>::decompose<float> (const SparseMatrix<float> &,
                                                   const double);
template void SparseMIC<double>::vmult<float> (Vector<float> &,
                                               const Vector<float> &) const;

template class SparseMIC<float>;
template void SparseMIC<float>::initialize<double> (const SparseMatrix<double> &,
						    const AdditionalData data);
template void SparseMIC<float>::decompose<double> (const SparseMatrix<double> &,
                                                   const double);
template void SparseMIC<float>::vmult<double> (Vector<double> &,
                                               const Vector<double> &) const;
template void SparseMIC<float>::initialize<float> (const SparseMatrix<float> &,
						   const AdditionalData data);
template void SparseMIC<float>::decompose<float> (const SparseMatrix<float> &,
                                                  const double);
template void SparseMIC<float>::vmult<float> (Vector<float> &,
                                              const Vector<float> &) const;


DEAL_II_NAMESPACE_CLOSE
