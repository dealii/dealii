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
#include <lac/full_matrix.h>
#include <lac/vector.h>

#include <iostream>

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
{
  state = LAPACKSupport::matrix;
}


template <typename number>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator = (const LAPACKFullMatrix<number>& M)
{
  TransposeTable<number>::operator=(M);
  state = LAPACKSupport::matrix;
  return *this;
}


template <typename number>
template <typename number2>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator = (const FullMatrix<number2>& M)
{
  Assert (this->n_rows() == M.m(), ExcDimensionMismatch(this->n_rows(), M.m()));
  Assert (this->n_cols() == M.n(), ExcDimensionMismatch(this->n_cols(), M.n()));
  for (unsigned int i=0;i<this->n_rows();++i)
    for (unsigned int j=0;j<this->n_cols();++j)
      (*this)(i,j) = M(i,j);
  
  state = LAPACKSupport::matrix;
  return *this;
}


template <typename number>
LAPACKFullMatrix<number> &
LAPACKFullMatrix<number>::operator = (const double d)
{
  Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());
  
  if (this->n_elements() != 0)
    std::fill_n (this->val, this->n_elements(), number());  
  
  state = LAPACKSupport::matrix;
  return *this;
}

#ifdef HAVE_LIBBLAS

template <typename number>
void
LAPACKFullMatrix<number>::vmult (
  Vector<number>       &w,
  const Vector<number> &v,
  const bool            adding) const
{
  Assert (state == matrix, ExcInvalidState());
  
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);
  
  gemv("N", &mm, &nn, &alpha, data(), &mm, v.val, &one, &beta, w.val, &one);
}


template <typename number>
void
LAPACKFullMatrix<number>::Tvmult (
  Vector<number>       &w,
  const Vector<number> &v,
  const bool            adding) const
{
  Assert (state == matrix, ExcInvalidState());
  
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  const number alpha = 1.;
  const number beta = (adding ? 1. : 0.);
  
  gemv("T", &mm, &nn, &alpha, data(), &mm, v.val, &one, &beta, w.val, &one);
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

#ifdef HAVE_LIBLAPACK

template <typename number>
void
LAPACKFullMatrix<number>::compute_eigenvalues()
{
  const int nn = this->n_cols();
  wr.resize(nn);
  wi.resize(nn);
  number* values = const_cast<number*> (this->data());
  
  int info;
  int lwork = -1;
  work.resize(1);
  geev(&N, &N, &nn, values, &nn,
       &wr[0], &wi[0],
       0, &one, 0, &one,
       &work[0], &lwork, &info);
				   // geev returns info=0 on
				   // success. Since we only queried
				   // the optimal size for work,
				   // everything else would not be
				   // acceptable.
  Assert (info == 0, ExcInternalError());

				   // Allocate working array according
				   // to suggestion.
  lwork = (int) (work[0]+.1);
  work.resize((unsigned int) lwork);
				   // Finally compute the eigenvalues.
  geev(&N, &N, &nn, values, &nn,
       &wr[0], &wi[0],
       0, &one, 0, &one,
       &work[0], &lwork, &info);
				   // Negative return value implies a
				   // wrong argument. This should be
				   // internal.
  Assert (info >=0, ExcInternalError());
//TODO:[GK] What if the QR method fails?
  if (info != 0)
    std::cerr << "LAPACK error in geev" << std::endl;

  state = LAPACKSupport::State(eigenvalues | unusable);
}

#endif

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


template class LAPACKFullMatrix<double>;
template LAPACKFullMatrix<double> &
LAPACKFullMatrix<double>::operator = (const FullMatrix<double>& M);
template LAPACKFullMatrix<double> &
LAPACKFullMatrix<double>::operator = (const FullMatrix<float>& M);
