//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/sparse_matrix.templates.h>
#include <lac/block_vector.h>

#define TYPEMAT double

template class SparseMatrix<TYPEMAT>;

#define TYPE2 float
#include "sparse_matrix_matrix.in.h"
#undef TYPE2

#define TYPE2 double
#include "sparse_matrix_matrix.in.h"
#undef TYPE2

#define TYPEVEC float
#include "sparse_matrix_vector.in.h"
#undef TYPEVEC

#define TYPEVEC double
#include "sparse_matrix_vector.in.h"
#undef TYPEVEC

				 // a prerelease of gcc3.0 fails to
				 // compile this due to long double
//  #undef TYPEVEC
//  #define TYPEVEC long double

//  #include <lac/sparse_matrix.2.templates>

#undef TYPEMAT


// instantiate a few more functions. this whole mess in this file (and the
// other sparse_matrix.* files) really needs to be cleaned up a little
template void SparseMatrix<double>::
vmult (Vector<double> &, const Vector<float> &) const;
template void SparseMatrix<double>::
Tvmult (Vector<double> &, const Vector<float> &) const;
template void SparseMatrix<double>::
vmult_add (Vector<double> &, const Vector<float> &) const;
template void SparseMatrix<double>::
Tvmult_add (Vector<double> &, const Vector<float> &) const;

template void SparseMatrix<double>::
vmult (Vector<float> &, const Vector<double> &) const;
template void SparseMatrix<double>::
Tvmult (Vector<float> &, const Vector<double> &) const;
template void SparseMatrix<double>::
vmult_add (Vector<float> &, const Vector<double> &) const;
template void SparseMatrix<double>::
Tvmult_add (Vector<float> &, const Vector<double> &) const;


template void SparseMatrix<double>::
vmult (BlockVector<double> &, const BlockVector<float> &) const;
template void SparseMatrix<double>::
Tvmult (BlockVector<double> &, const BlockVector<float> &) const;
template void SparseMatrix<double>::
vmult_add (BlockVector<double> &, const BlockVector<float> &) const;
template void SparseMatrix<double>::
Tvmult_add (BlockVector<double> &, const BlockVector<float> &) const;

template void SparseMatrix<double>::
vmult (BlockVector<float> &, const BlockVector<double> &) const;
template void SparseMatrix<double>::
Tvmult (BlockVector<float> &, const BlockVector<double> &) const;
template void SparseMatrix<double>::
vmult_add (BlockVector<float> &, const BlockVector<double> &) const;
template void SparseMatrix<double>::
Tvmult_add (BlockVector<float> &, const BlockVector<double> &) const;


template void SparseMatrix<double>::
vmult (Vector<double> &, const BlockVector<float> &) const;
template void SparseMatrix<double>::
Tvmult (Vector<double> &, const BlockVector<float> &) const;
template void SparseMatrix<double>::
vmult_add (Vector<double> &, const BlockVector<float> &) const;
template void SparseMatrix<double>::
Tvmult_add (Vector<double> &, const BlockVector<float> &) const;

template void SparseMatrix<double>::
vmult (BlockVector<float> &, const Vector<double> &) const;
template void SparseMatrix<double>::
Tvmult (BlockVector<float> &, const Vector<double> &) const;
template void SparseMatrix<double>::
vmult_add (BlockVector<float> &, const Vector<double> &) const;
template void SparseMatrix<double>::
Tvmult_add (BlockVector<float> &, const Vector<double> &) const;


template void SparseMatrix<double>::
vmult (BlockVector<double> &, const Vector<float> &) const;
template void SparseMatrix<double>::
Tvmult (BlockVector<double> &, const Vector<float> &) const;
template void SparseMatrix<double>::
vmult_add (BlockVector<double> &, const Vector<float> &) const;
template void SparseMatrix<double>::
Tvmult_add (BlockVector<double> &, const Vector<float> &) const;

template void SparseMatrix<double>::
vmult (Vector<float> &, const BlockVector<double> &) const;
template void SparseMatrix<double>::
Tvmult (Vector<float> &, const BlockVector<double> &) const;
template void SparseMatrix<double>::
vmult_add (Vector<float> &, const BlockVector<double> &) const;
template void SparseMatrix<double>::
Tvmult_add (Vector<float> &, const BlockVector<double> &) const;


template void SparseMatrix<double>::
vmult (Vector<double> &, const BlockVector<double> &) const;
template void SparseMatrix<double>::
Tvmult (Vector<double> &, const BlockVector<double> &) const;
template void SparseMatrix<double>::
vmult_add (Vector<double> &, const BlockVector<double> &) const;
template void SparseMatrix<double>::
Tvmult_add (Vector<double> &, const BlockVector<double> &) const;

template void SparseMatrix<double>::
vmult (BlockVector<double> &, const Vector<double> &) const;
template void SparseMatrix<double>::
Tvmult (BlockVector<double> &, const Vector<double> &) const;
template void SparseMatrix<double>::
vmult_add (BlockVector<double> &, const Vector<double> &) const;
template void SparseMatrix<double>::
Tvmult_add (BlockVector<double> &, const Vector<double> &) const;
