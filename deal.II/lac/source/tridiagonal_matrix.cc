//-------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998 - 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------

#include <lac/tridiagonal_matrix.h>
#include <lac/vector.h>

template<typename number>
TridiagonalMatrix<number>::TridiagonalMatrix(unsigned int size)
		:
		diagonal(size, 0.),
		left(size, 0.),
		right(size, 0.)
{}


template<typename number>
bool
TridiagonalMatrix<number>::all_zero() const
{
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
  Assert(w.size() == n(), ExcDimensionMismatch(w.size(), n()));
  Assert(v.size() == n(), ExcDimensionMismatch(v.size(), n()));

  if (n()==0) return;
  
  if (adding)
    {
      w(0) += diagonal[0]*v(0) + right[0]*v(1);
      const unsigned int e=n()-1;
      for (unsigned int i=1;i<e;++i)
	w(i) += left[i]*v(i-1)+diagonal[i]*v(i)+right[i]*v(i+1);
      w(e) += left[e]*v(e-1)+diagonal[e]*v(e);
    }
  else
    {
      w(0) = diagonal[0]*v(0) + right[0]*v(1);
      const unsigned int e=n()-1;
      for (unsigned int i=1;i<e;++i)
	w(i) = left[i]*v(i-1)+diagonal[i]*v(i)+right[i]*v(i+1);
      w(e) = left[e]*v(e-1)+diagonal[e]*v(e);
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
  Assert(w.size() == n(), ExcDimensionMismatch(w.size(), n()));
  Assert(v.size() == n(), ExcDimensionMismatch(v.size(), n()));

  if (n()==0) return;
//TODO:[GK] Check this!!!  
  if (adding)
    {
      w(0) += diagonal[0]*v(0) + left[1]*v(1);
      const unsigned int e=n()-1;
      for (unsigned int i=1;i<e;++i)
	w(i) += left[i+1]*v(i+1)+diagonal[i]*v(i)+right[i-1]*v(i-1);
      w(e) += right[e-1]*v(e-1)+diagonal[e]*v(e);
    }
  else
    {
      w(0) = diagonal[0]*v(0) + left[1]*v(1);
      const unsigned int e=n()-1;
      for (unsigned int i=1;i<e;++i)
	w(i) = left[i+1]*v(i+1)+diagonal[i]*v(i)+right[i-1]*v(i-1);
      w(e) = right[e-1]*v(e-1)+diagonal[e]*v(e);
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

/*
template<typename number>
TridiagonalMatrix<number>::
{
}


template<typename number>
TridiagonalMatrix<number>::
{
}


*/

template TridiagonalMatrix<float>;
template TridiagonalMatrix<double>;
