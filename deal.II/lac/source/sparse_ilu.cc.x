//----------------------------  sparse_ilu.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_ilu.cc  ---------------------------


#include <lac/sparse_ilu.templates.h>


// explicit instantiations
template class SparseILU<double>;
template void SparseILU<double>::decompose<double> (const SparseMatrix<double> &,
						    const double);
template void SparseILU<double>::apply_decomposition<double> (Vector<double> &,
							      const Vector<double> &) const;
template void SparseILU<double>::decompose<float> (const SparseMatrix<float> &,
						   const double);
template void SparseILU<double>::apply_decomposition<float> (Vector<float> &,
							     const Vector<float> &) const;


template class SparseILU<float>;
template void SparseILU<float>::decompose<double> (const SparseMatrix<double> &,
						   const double);
template void SparseILU<float>::apply_decomposition<double> (Vector<double> &,
							     const Vector<double> &) const;
template void SparseILU<float>::decompose<float> (const SparseMatrix<float> &,
						  const double);
template void SparseILU<float>::apply_decomposition<float> (Vector<float> &,
							    const Vector<float> &) const;
