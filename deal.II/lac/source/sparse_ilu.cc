//---------------------------------------------------------------------------
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/sparse_ilu.templates.h>

DEAL_II_NAMESPACE_OPEN


// explicit instantiations
template class SparseILU<double>;
template void SparseILU<double>::initialize<double> (const SparseMatrix<double> &,
						     const AdditionalData data);
template void SparseILU<double>::decompose<double> (const SparseMatrix<double> &,
						    const double);
template void SparseILU<double>::vmult <double> (Vector<double> &,
                                                 const Vector<double> &) const;
template void SparseILU<double>::Tvmult <double> (Vector<double> &,
						  const Vector<double> &) const;
template void SparseILU<double>::initialize<float> (const SparseMatrix<float> &,
						    const AdditionalData data);
template void SparseILU<double>::decompose<float> (const SparseMatrix<float> &,
						   const double);
template void SparseILU<double>::vmult<float> (Vector<float> &,
                                               const Vector<float> &) const;
template void SparseILU<double>::Tvmult<float> (Vector<float> &,
						const Vector<float> &) const;


template class SparseILU<float>;
template void SparseILU<float>::initialize<double> (const SparseMatrix<double> &,
						    const AdditionalData data);
template void SparseILU<float>::decompose<double> (const SparseMatrix<double> &,
						   const double);
template void SparseILU<float>::vmult<double> (Vector<double> &,
                                               const Vector<double> &) const;
template void SparseILU<float>::Tvmult<double> (Vector<double> &,
						const Vector<double> &) const;
template void SparseILU<float>::initialize<float> (const SparseMatrix<float> &,
						   const AdditionalData data);
template void SparseILU<float>::decompose<float> (const SparseMatrix<float> &,
						  const double);
template void SparseILU<float>::vmult<float> (Vector<float> &,
                                              const Vector<float> &) const;
template void SparseILU<float>::Tvmult<float> (Vector<float> &,
					       const Vector<float> &) const;

DEAL_II_NAMESPACE_CLOSE
