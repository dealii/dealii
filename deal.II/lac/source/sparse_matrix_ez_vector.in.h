//----------------------------  sparse_matrix_ez.2.templates  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_matrix_ez.2.templates  ---------------------------


// File copied from sparse_matrix.2.templates and modified

// Driver file for SparseMatrixEZ functions with two types.

// TYPEMAT and TYPEVEC are defined in sparse_matrix_ez?.cc

template void SparseMatrixEZ<TYPEMAT>::vmult<TYPEVEC> (
  Vector<TYPEVEC> &,
  const Vector<TYPEVEC> &) const;
template void SparseMatrixEZ<TYPEMAT>::Tvmult<TYPEVEC> (
  Vector<TYPEVEC> &,
  const Vector<TYPEVEC> &) const;
template void SparseMatrixEZ<TYPEMAT>::vmult_add<TYPEVEC> (
  Vector<TYPEVEC> &,
  const Vector<TYPEVEC> &) const;
template void SparseMatrixEZ<TYPEMAT>::Tvmult_add<TYPEVEC> (
  Vector<TYPEVEC> &,
  const Vector<TYPEVEC> &) const;

//template TYPEVEC
//SparseMatrixEZ<TYPEMAT>::matrix_norm_square<TYPEVEC> (const Vector<TYPEVEC> &) const;

//template TYPEVEC
//SparseMatrixEZ<TYPEMAT>::matrix_scalar_product<TYPEVEC> (const Vector<TYPEVEC> &,
//						       const Vector<TYPEVEC> &) const;

//template TYPEVEC SparseMatrixEZ<TYPEMAT>::residual<TYPEVEC> (Vector<TYPEVEC> &,
//							 const Vector<TYPEVEC> &,
//							 const Vector<TYPEVEC> &) const;

template void SparseMatrixEZ<TYPEMAT>::precondition_SSOR<TYPEVEC> (
  Vector<TYPEVEC> &,
  const Vector<TYPEVEC> &,
  const TYPEMAT) const;

template void SparseMatrixEZ<TYPEMAT>::precondition_SOR<TYPEVEC> (
  Vector<TYPEVEC> &,
  const Vector<TYPEVEC> &,
  const TYPEMAT) const;

template void SparseMatrixEZ<TYPEMAT>::precondition_TSOR<TYPEVEC> (
  Vector<TYPEVEC> &,
  const Vector<TYPEVEC> &,
  const TYPEMAT) const;

template void SparseMatrixEZ<TYPEMAT>::precondition_Jacobi<TYPEVEC> (
  Vector<TYPEVEC> &,
  const Vector<TYPEVEC> &,
  const TYPEMAT) const;

