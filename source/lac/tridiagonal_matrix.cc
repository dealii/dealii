// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


#include <deal.II/lac/tridiagonal_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/lapack_templates.h>

DEAL_II_NAMESPACE_OPEN

using namespace LAPACKSupport;

template<typename number>
TridiagonalMatrix<number>::TridiagonalMatrix(
  size_type size,
  bool symmetric)
  :
  diagonal(size, 0.),
  left((symmetric ? 0 : size), 0.),
  right(size, 0.),
  is_symmetric(symmetric),
  state(matrix)
{}


template<typename number>
void
TridiagonalMatrix<number>::reinit(
  size_type size,
  bool symmetric)
{
  is_symmetric = symmetric;
  diagonal.resize(size);
  right.resize(size);
  left.resize(symmetric ? 0 : size);
  state = matrix;
}


template<typename number>
bool
TridiagonalMatrix<number>::all_zero() const
{
  Assert(state == matrix, ExcState(state));

  typename std::vector<number>::const_iterator i;
  typename std::vector<number>::const_iterator e;

  e = diagonal.end();
  for (i=diagonal.begin() ; i != e ; ++i)
    if (*i != 0.) return false;

  e = left.end();
  for (i=left.begin() ; i != e ; ++i)
    if (*i != 0.) return false;

  e = right.end();
  for (i=right.begin() ; i != e ; ++i)
    if (*i != 0.) return false;
  return true;
}


template<typename number>
void
TridiagonalMatrix<number>::vmult (
  Vector<number>       &w,
  const Vector<number> &v,
  const bool            adding) const
{
  Assert(state == matrix, ExcState(state));

  Assert(w.size() == n(), ExcDimensionMismatch(w.size(), n()));
  Assert(v.size() == n(), ExcDimensionMismatch(v.size(), n()));

  if (n()==0) return;

  // The actual loop skips the first
  // and last row
  const size_type e=n()-1;
  // Let iterators point to the first
  // entry of each diagonal
  typename std::vector<number>::const_iterator d = diagonal.begin();
  typename std::vector<number>::const_iterator r = right.begin();
  // The left diagonal starts one
  // later or is equal to the right
  // one for symmetric storage
  typename std::vector<number>::const_iterator l = left.begin();
  if (is_symmetric)
    l = r;
  else
    ++l;

  if (adding)
    {
      // Treat first row separately
      w(0) += (*d) * v(0) + (*r) * v(1);
      ++d;
      ++r;
      // All rows with three entries
      for (size_type i=1; i<e; ++i,++d,++r,++l)
        w(i) += (*l) * v(i-1) + (*d) * v(i) + (*r) * v(i+1);
      // Last row is special again
      w(e) += (*l) * v(e-1) + (*d) * v(e);
    }
  else
    {
      w(0) = (*d) * v(0) + (*r) * v(1);
      ++d;
      ++r;
      for (size_type i=1; i<e; ++i,++d,++r,++l)
        w(i) = (*l) * v(i-1) + (*d) * v(i) + (*r) * v(i+1);
      w(e) = (*l) * v(e-1) + (*d) * v(e);
    }
}


template<typename number>
void
TridiagonalMatrix<number>::vmult_add (
  Vector<number>       &w,
  const Vector<number> &v) const
{
  vmult(w, v, true);
}


template<typename number>
void
TridiagonalMatrix<number>::Tvmult (
  Vector<number>       &w,
  const Vector<number> &v,
  const bool            adding) const
{
  Assert(state == matrix, ExcState(state));

  Assert(w.size() == n(), ExcDimensionMismatch(w.size(), n()));
  Assert(v.size() == n(), ExcDimensionMismatch(v.size(), n()));

  if (n()==0) return;

  const size_type e=n()-1;
  typename std::vector<number>::const_iterator d = diagonal.begin();
  typename std::vector<number>::const_iterator r = right.begin();
  typename std::vector<number>::const_iterator l = left.begin();
  if (is_symmetric)
    l = r;
  else
    ++l;

  if (adding)
    {
      w(0) += (*d) * v(0) + (*l) * v(1);
      ++d;
      ++l;
      for (size_type i=1; i<e; ++i,++d,++r,++l)
        w(i) += (*l) * v(i+1) + (*d) * v(i) + (*r) * v(i-1);
      w(e) += (*d) * v(e) + (*r) * v(e-1);
    }
  else
    {
      w(0) = (*d) * v(0) + (*l) * v(1);
      ++d;
      ++l;
      for (size_type i=1; i<e; ++i,++d,++r,++l)
        w(i) = (*l) * v(i+1) + (*d) * v(i) + (*r) * v(i-1);
      w(e) = (*d) * v(e) + (*r) * v(e-1);
    }
}


template<typename number>
void
TridiagonalMatrix<number>::Tvmult_add (
  Vector<number>       &w,
  const Vector<number> &v) const
{
  Tvmult(w, v, true);
}


template<typename number>
number
TridiagonalMatrix<number>::matrix_scalar_product(
  const Vector<number> &w,
  const Vector<number> &v) const
{
  Assert(state == matrix, ExcState(state));

  const size_type e=n()-1;
  typename std::vector<number>::const_iterator d = diagonal.begin();
  typename std::vector<number>::const_iterator r = right.begin();
  typename std::vector<number>::const_iterator l = left.begin();
  if (is_symmetric)
    l = r;
  else
    ++l;

  number result = w(0) * ((*d) * v(0) + (*r) * v(1));
  ++d;
  ++r;
  for (size_type i=1; i<e; ++i,++d,++r,++l)
    result += w(i) * ((*l) * v(i-1)+ (*d) * v(i)+ (*r) * v(i+1));
  result += w(e) * ((*l) * v(e-1) + (*d) * v(e));
  return result;
}


template<typename number>
number
TridiagonalMatrix<number>::matrix_norm_square(
  const Vector<number> &v) const
{
  return matrix_scalar_product(v,v);
}


template<>
void
TridiagonalMatrix<double>::compute_eigenvalues()
{
#ifdef DEAL_II_WITH_LAPACK
  Assert(state == matrix, ExcState(state));
  Assert(is_symmetric, ExcNotImplemented());

  const int nn = n();
  int info;
  stev (&N, &nn, &*diagonal.begin(), &*right.begin(), 0, &one, 0, &info);
  Assert(info == 0, ExcInternalError());

  state = eigenvalues;
#else
  Assert(false, ExcNeedsLAPACK());
#endif
}


template<typename number>
number
TridiagonalMatrix<number>::eigenvalue(const size_type i) const
{
  Assert(state == eigenvalues, ExcState(state));
  Assert(i<n(), ExcIndexRange(i,0,n()));
  return diagonal[i];
}


/*
template<typename number>
TridiagonalMatrix<number>::
{
}


*/

template class TridiagonalMatrix<float>;
template class TridiagonalMatrix<double>;

DEAL_II_NAMESPACE_CLOSE
