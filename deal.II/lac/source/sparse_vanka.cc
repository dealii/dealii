//----------------------------  sparse_vanka.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_vanka.cc  ---------------------------


#include <lac/sparse_vanka.templates.h>


// explicit instantiations
template class SparseVanka<float>;
template class SparseVanka<double>;

template void SparseVanka<double>::operator () (Vector<float>       &dst,
						const Vector<float> &src) const;
template void SparseVanka<double>::operator () (Vector<double>       &dst,
						const Vector<double> &src) const;


template class SparseBlockVanka<float>;
template class SparseBlockVanka<double>;

template void SparseBlockVanka<double>::operator () (Vector<float>       &dst,
						     const Vector<float> &src) const;
template void SparseBlockVanka<double>::operator () (Vector<double>       &dst,
						     const Vector<double> &src) const;
