//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/lapack_full_matrix.h>
#include <lac/lapack_templates.h>
#include <lac/lapack_support.h>
#include <lac/vector.h>

using namespace LAPACKSupport;

template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix(const unsigned int n)
		:
		TransposeTable<number> (n,n)
{}


template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix(
  const unsigned int m,
  const unsigned int n)
		:
		TransposeTable<number> (m,n)
{}


template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix(const LAPACKFullMatrix &M)
		:
		TransposeTable<number> (M)
{}


template <typename number>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator = (const LAPACKFullMatrix<number>& M)
{
  TransposeTable<number>::operator=(M);
  return *this;
}


template <typename number>
template <typename number2>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator = (const FullMatrix<number2>& M)
{
  Assert (n() == M.m(), ExcDimensionMismatch(n(), M.m()));
  Assert (n() == M.n(), ExcDimensionMismatch(n(), M.n()));
  for (unsigned int i=0;i<m;++i)
    for (unsigned int j=0;j<n;++j)
      (*this)(i,j) = M(i,j);
  return *this;
}


template <typename number>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator = (const double d)
{
  Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());
  
  if (this->n_elements() != 0)
    std::fill_n (this->val, this->n_elements(), number());
  
  return *this;
}

#ifdef HAVE_LIBLAPACK


template <typename number>
void
LAPACKFullMatrix<number>::vmult (
  Vector<number>       &w,
  const Vector<number> &v,
  const bool            adding) const
{
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);
  
  gemv(&N, &mm, &nn, &alpha, data(), &mm, v.val, &one, &beta, w.val, &one);
}


template <typename number>
void
LAPACKFullMatrix<number>::Tvmult (
  Vector<number>       &w,
  const Vector<number> &v,
  const bool            adding) const
{
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);
  
  gemv(&T, &mm, &nn, &alpha, data(), &mm, v.val, &one, &beta, w.val, &one);
}

#else


template <typename number>
void
LAPACKFullMatrix<number>::vmult (
  Vector<number>       &,
  const Vector<number> &,
  const bool            ) const
{
  Assert(false, ExcNeedsBLAS());
}


template <typename number>
void
LAPACKFullMatrix<number>::Tvmult (
  Vector<number>       &,
  const Vector<number> &,
  const bool            ) const
{
  Assert(false, ExcNeedsBLAS());
}


#endif

// template <typename number>
// LAPACKFullMatrix<number>::()
// {}


// template <typename number>
// LAPACKFullMatrix<number>::()
// {}


// template <typename number>
// LAPACKFullMatrix<number>::()
// {}


template <typename number>
void
LAPACKFullMatrix<number>::vmult_add (
  Vector<number>       &w,
  const Vector<number> &v) const
{
  vmult(w, v, true);
}


template <typename number>
void
LAPACKFullMatrix<number>::Tvmult_add (
  Vector<number>       &w,
  const Vector<number> &v) const
{
  Tvmult(w, v, true);
}


template LAPACKFullMatrix<double>;
