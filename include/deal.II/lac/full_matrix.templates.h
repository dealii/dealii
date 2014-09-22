// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__full_matrix_templates_h
#define __deal2__full_matrix_templates_h


//TODO: this file has a lot of operations between matrices and matrices or matrices and vectors of different precision. we should go through the file and in each case pick the more accurate data type for intermediate results. currently, the choice is pretty much random. this may also allow us some operations where one operand is complex and the other is not

#include <deal.II/base/config.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>

#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN


template <typename number>
FullMatrix<number>::FullMatrix (const size_type n)
  :
  Table<2,number> (n,n)
{}


template <typename number>
FullMatrix<number>::FullMatrix (const size_type m,
                                const size_type n)
  :
  Table<2,number> (m, n)
{}


template <typename number>
FullMatrix<number>::FullMatrix (const size_type m,
                                const size_type n,
                                const number *entries)
  :
  Table<2,number> (m, n)
{
  this->fill (entries);
}


template <typename number>
FullMatrix<number>::FullMatrix (const FullMatrix &m)
  :
  Table<2,number> (m)
{}



template <typename number>
FullMatrix<number>::FullMatrix (const IdentityMatrix &id)
  :
  Table<2,number> (id.m(), id.n())
{
  for (size_type i=0; i<id.m(); ++i)
    (*this)(i,i) = 1;
}


template <typename number>
FullMatrix<number> &
FullMatrix<number>::operator = (const FullMatrix<number> &M)
{
  Table<2,number>::operator=(M);
  return *this;
}


template <typename number>
template <typename number2>
FullMatrix<number> &
FullMatrix<number>::operator = (const FullMatrix<number2> &M)
{
  TableBase<2,number>::operator=(M);
  return *this;
}



template <typename number>
FullMatrix<number> &
FullMatrix<number>::operator = (const IdentityMatrix &id)
{
  this->reinit (id.m(), id.n());
  for (size_type i=0; i<id.m(); ++i)
    (*this)(i,i) = 1.;

  return *this;
}



template <typename number>
template <typename number2>
FullMatrix<number> &
FullMatrix<number>::operator = (const LAPACKFullMatrix<number2> &M)
{
  Assert (this->m() == M.n_rows(), ExcDimensionMismatch(this->m(), M.n_rows()));
  Assert (this->n() == M.n_cols(), ExcDimensionMismatch(this->n(), M.n_rows()));
  for (size_type i=0; i<this->m(); ++i)
    for (size_type j=0; j<this->n(); ++j)
      (*this)(i,j) = M(i,j);

  return *this;
}



template <typename number>
bool
FullMatrix<number>::all_zero () const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  const number *p = &this->values[0];
  const number *const e = &this->values[0] + this->n_elements();
  while (p!=e)
    if (*p++ != number(0.0))
      return false;

  return true;
}



template <typename number>
FullMatrix<number> &
FullMatrix<number>::operator *= (const number factor)
{

  Assert (numbers::is_finite(factor), ExcNumberNotFinite());

  number       *p = &(*this)(0,0);
  const number *e = &(*this)(0,0) + n()*m();
  while (p != e)
    *p++ *= factor;

  return *this;
}



template <typename number>
FullMatrix<number> &
FullMatrix<number>::operator /= (const number factor)
{

  Assert (numbers::is_finite(factor), ExcNumberNotFinite());

  number       *p = &(*this)(0,0);
  const number *e = &(*this)(0,0) + n()*m();

  const number factor_inv = number(1.)/factor;

  Assert (numbers::is_finite(factor_inv), ExcNumberNotFinite());

  while (p != e)
    *p++ *= factor_inv;

  return *this;
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::vmult (Vector<number2> &dst,
                           const Vector<number2> &src,
                           const bool adding) const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert(dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));

  Assert (&src != &dst, ExcSourceEqualsDestination());

  const number *e = &this->values[0];
  // get access to the data in order to
  // avoid copying it when using the ()
  // operator
  const number2 *src_ptr = &(*const_cast<Vector<number2>*>(&src))(0);
  const size_type size_m = m(), size_n = n();
  for (size_type i=0; i<size_m; ++i)
    {
      number2 s = adding ? dst(i) : 0.;
      for (size_type j=0; j<size_n; ++j)
        s += src_ptr[j] * number2(*(e++));
      dst(i) = s;
    }
}



template <typename number>
template <typename number2>
void FullMatrix<number>::Tvmult (Vector<number2>       &dst,
                                 const Vector<number2> &src,
                                 const bool             adding) const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert(dst.size() == n(), ExcDimensionMismatch(dst.size(), n()));
  Assert(src.size() == m(), ExcDimensionMismatch(src.size(), m()));

  Assert (&src != &dst, ExcSourceEqualsDestination());

  const number *e = &this->values[0];
  number2 *dst_ptr = &dst(0);
  const size_type size_m = m(), size_n = n();

  // zero out data if we are not adding
  if (!adding)
    for (size_type j=0; j<size_n; ++j)
      dst_ptr[j] = 0.;

  // write the loop in a way that we can
  // access the data contiguously
  for (size_type i=0; i<size_m; ++i)
    {
      const number2 d = src(i);
      for (size_type j=0; j<size_n; ++j)
        dst_ptr[j] += d * number2(*(e++));
    };
}


template <typename number>
template <typename number2, typename number3>
number FullMatrix<number>::residual (Vector<number2> &dst,
                                     const Vector<number2> &src,
                                     const Vector<number3> &right) const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert(dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));
  Assert(right.size() == m(), ExcDimensionMismatch(right.size(), m()));

  Assert (&src != &dst, ExcSourceEqualsDestination());

  number res = 0.;
  const size_type size_m = m(),
                  size_n = n();
  for (size_type i=0; i<size_n; ++i)
    {
      number s = number(right(i));
      for (size_type j=0; j<size_m; ++j)
        s -= number(src(j)) * (*this)(i,j);
      dst(i) = s;
      res += s*s;
    }
  return std::sqrt(res);
}



template <typename number>
template <typename number2>
void FullMatrix<number>::forward (Vector<number2>       &dst,
                                  const Vector<number2> &src) const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));
  Assert (src.size() == n(), ExcDimensionMismatch(src.size(), n()));

  size_type i,j;
  size_type nu = ( (m()<n()) ? m() : n());
  for (i=0; i<nu; ++i)
    {
      number s = number(src(i));
      for (j=0; j<i; ++j)
        s -= number(dst(j)) * (*this)(i,j);
      dst(i) = s/(*this)(i,i);
      Assert(numbers::is_finite(dst(i)), ExcNumberNotFinite());
    }
}



template <typename number>
template <typename number2>
void FullMatrix<number>::backward (Vector<number2>       &dst,
                                   const Vector<number2> &src) const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  size_type j;
  size_type nu = (m()<n() ? m() : n());
  for (int i=nu-1; i>=0; --i)
    {
      number2 s = src(i);
      for (j=i+1; j<nu; ++j)
        s -= dst(j) * number2((*this)(i,j));
      dst(i) = s/number2((*this)(i,i));
      Assert(numbers::is_finite(dst(i)), ExcNumberNotFinite());
    }
}



template <typename number>
template <typename number2>
void FullMatrix<number>::fill (const FullMatrix<number2> &src,
                               const size_type dst_offset_i,
                               const size_type dst_offset_j,
                               const size_type src_offset_i,
                               const size_type src_offset_j)
{
  Assert (dst_offset_i < m(),
          ExcIndexRange (dst_offset_i, 0, m()));
  Assert (dst_offset_j < n(),
          ExcIndexRange (dst_offset_j, 0, n()));
  Assert (src_offset_i < src.m(),
          ExcIndexRange (src_offset_i, 0, src.m()));
  Assert (src_offset_j < src.n(),
          ExcIndexRange (src_offset_j, 0, src.n()));

  // Compute maximal size of copied block
  const size_type rows = std::min (m() - dst_offset_i,
                                   src.m() - src_offset_i);
  const size_type cols = std::min (n() - dst_offset_j,
                                   src.n() - src_offset_j);

  for (size_type i=0; i<rows ; ++i)
    for (size_type j=0; j<cols ; ++j)
      (*this)(dst_offset_i+i,dst_offset_j+j)
        = src(src_offset_i+i,src_offset_j+j);
}


template <typename number>
template <typename number2>
void FullMatrix<number>::fill_permutation (const FullMatrix<number2> &src,
                                           const std::vector<size_type> &p_rows,
                                           const std::vector<size_type> &p_cols)
{
  Assert (p_rows.size() == this->n_rows(),
          ExcDimensionMismatch (p_rows.size(), this->n_rows()));
  Assert (p_cols.size() == this->n_cols(),
          ExcDimensionMismatch (p_cols.size(), this->n_cols()));

  for (size_type i=0; i<this->n_rows(); ++i)
    for (size_type j=0; j<this->n_cols(); ++j)
      (*this)(i,j) = src(p_rows[i], p_cols[j]);
}



/*  template <typename number> */
/*  template <typename number2> */
/*  void FullMatrix<number>::fill (const number2* entries) */
/*  { */
/*      if (n_cols()*n_rows() != 0) */
/*        std::copy (entries, entries+n_rows()*n_cols(), &this->values[0]); */
/*  } */



template <typename number>
void FullMatrix<number>::add_row (const size_type i,
                                  const number s,
                                  const size_type j)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  for (size_type k=0; k<m(); ++k)
    (*this)(i,k) += s*(*this)(j,k);
}


template <typename number>
void FullMatrix<number>::add_row (const size_type i,
                                  const number s,
                                  const size_type j,
                                  const number t,
                                  const size_type k)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  const size_type size_m = m();
  for (size_type l=0; l<size_m; ++l)
    (*this)(i,l) += s*(*this)(j,l) + t*(*this)(k,l);
}


template <typename number>
void FullMatrix<number>::add_col (const size_type i, const number s,
                                  const size_type j)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  for (size_type k=0; k<n(); ++k)
    (*this)(k,i) += s*(*this)(k,j);
}


template <typename number>
void FullMatrix<number>::add_col (const size_type i, const number s,
                                  const size_type j, const number t,
                                  const size_type k)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  for (size_t l=0; l<n(); ++l)
    (*this)(l,i) += s*(*this)(l,j) + t*(*this)(l,k);
}



template <typename number>
void FullMatrix<number>::swap_row (const size_type i,
                                   const size_type j)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  for (size_type k=0; k<n(); ++k)
    std::swap ((*this)(i,k),
               (*this)(j,k));
}


template <typename number>
void FullMatrix<number>::swap_col (const size_type i,
                                   const size_type j)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  for (size_type k=0; k<m(); ++k)
    std::swap ((*this)(k,i),
               (*this)(k,j));
}


template <typename number>
void FullMatrix<number>::diagadd (const number src)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  Assert (m() == n(), ExcDimensionMismatch(m(),n()));

  for (size_type i=0; i<n(); ++i)
    (*this)(i,i) += src;
}


template <typename number>
template <typename number2>
void FullMatrix<number>::equ (const number               a,
                              const FullMatrix<number2> &A)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (m() == A.m(), ExcDimensionMismatch(m(), A.m()));
  Assert (n() == A.n(), ExcDimensionMismatch(n(), A.n()));

  for (size_type i=0; i<m(); ++i)
    for (size_type j=0; j<n(); ++j)
      (*this)(i,j) = a * number(A(i,j));
}


template <typename number>
template <typename number2>
void
FullMatrix<number>::equ (const number               a,
                         const FullMatrix<number2> &A,
                         const number               b,
                         const FullMatrix<number2> &B)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (m() == A.m(), ExcDimensionMismatch(m(), A.m()));
  Assert (n() == A.n(), ExcDimensionMismatch(n(), A.n()));
  Assert (m() == B.m(), ExcDimensionMismatch(m(), B.m()));
  Assert (n() == B.n(), ExcDimensionMismatch(n(), B.n()));

  for (size_type i=0; i<m(); ++i)
    for (size_type j=0; j<n(); ++j)
      (*this)(i,j) = a * number(A(i,j)) + b * number(B(i,j));
}


template <typename number>
template <typename number2>
void
FullMatrix<number>::equ (const number               a,
                         const FullMatrix<number2> &A,
                         const number               b,
                         const FullMatrix<number2> &B,
                         const number               c,
                         const FullMatrix<number2> &C)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (m() == A.m(), ExcDimensionMismatch(m(), A.m()));
  Assert (n() == A.n(), ExcDimensionMismatch(n(), A.n()));
  Assert (m() == B.m(), ExcDimensionMismatch(m(), B.m()));
  Assert (n() == B.n(), ExcDimensionMismatch(n(), B.n()));
  Assert (m() == C.m(), ExcDimensionMismatch(m(), C.m()));
  Assert (n() == C.n(), ExcDimensionMismatch(n(), C.n()));

  for (size_type i=0; i<m(); ++i)
    for (size_type j=0; j<n(); ++j)
      (*this)(i,j) = a * number(A(i,j)) +
                     b * number(B(i,j)) +
                     c * number(C(i,j));
}



template <typename number>
template <typename number2>
void FullMatrix<number>::mmult (FullMatrix<number2>       &dst,
                                const FullMatrix<number2> &src,
                                const bool                 adding) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  Assert (n() == src.m(), ExcDimensionMismatch(n(), src.m()));
  Assert (dst.n() == src.n(), ExcDimensionMismatch(dst.n(), src.n()));
  Assert (dst.m() == m(), ExcDimensionMismatch(m(), dst.m()));

  // see if we can use BLAS algorithms for this and if the type for 'number'
  // works for us (it is usually not efficient to use BLAS for very small
  // matrices):
#if defined(HAVE_DGEMM_) && defined (HAVE_SGEMM_)
  if ((types_are_equal<number,double>::value
       ||
       types_are_equal<number,float>::value)
      &&
      types_are_equal<number,number2>::value)
    if (this->n()*this->m()*src.n() > 300)
      {
        // In case we have the BLAS function gemm detected at configure, we
        // use that algorithm for matrix-matrix multiplication since it
        // provides better performance than the deal.II native function (it
        // uses cache and register blocking in order to access local data).
        //
        // Note that BLAS/LAPACK stores matrix elements column-wise (i.e., all
        // values in one column, then all in the next, etc.), whereas the
        // FullMatrix stores them row-wise.  We ignore that difference, and
        // give our row-wise data to BLAS, let BLAS build the product of
        // transpose matrices, and read the result as if it were row-wise
        // again. In other words, we calculate (B^T A^T)^T, which is AB.

        const int m = src.n();
        const int n = this->m();
        const int k = this->n();
        const char *notrans = "n";

        const number alpha = 1.;
        const number beta = (adding == true) ? 1. : 0.;

        // Use the BLAS function gemm for calculating the matrix-matrix
        // product.
        gemm(notrans, notrans, &m, &n, &k, &alpha, &src(0,0), &m,
             &this->values[0], &k, &beta, &dst(0,0), &m);

        return;
      }

#endif

  const size_type m = this->m(), n = src.n(), l = this->n();

  // arrange the loops in a way that we keep write operations low, (writing is
  // usually more costly than reading), even though we need to access the data
  // in src not in a contiguous way.
  for (size_type i=0; i<m; i++)
    for (size_type j=0; j<n; j++)
      {
        number2 add_value = adding ? dst(i,j) : 0.;
        for (size_type k=0; k<l; k++)
          add_value += (number2)(*this)(i,k) * (number2)(src(k,j));
        dst(i,j) = add_value;
      }
}



template <typename number>
template <typename number2>
void FullMatrix<number>::Tmmult (FullMatrix<number2>       &dst,
                                 const FullMatrix<number2> &src,
                                 const bool                 adding) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == dst.m(), ExcDimensionMismatch(n(), dst.m()));
  Assert (src.n() == dst.n(), ExcDimensionMismatch(src.n(), dst.n()));


  // see if we can use BLAS algorithms for this and if the type for 'number'
  // works for us (it is usually not efficient to use BLAS for very small
  // matrices):
#if defined(HAVE_DGEMM_) && defined (HAVE_SGEMM_)
  if ((types_are_equal<number,double>::value
       ||
       types_are_equal<number,float>::value)
      &&
      types_are_equal<number,number2>::value)
    if (this->n()*this->m()*src.n() > 300)
      {
        // In case we have the BLAS function gemm detected at configure, we
        // use that algorithm for matrix-matrix multiplication since it
        // provides better performance than the deal.II native function (it
        // uses cache and register blocking in order to access local data).
        //
        // Note that BLAS/LAPACK stores matrix elements column-wise (i.e., all
        // values in one column, then all in the next, etc.), whereas the
        // FullMatrix stores them row-wise.  We ignore that difference, and
        // give our row-wise data to BLAS, let BLAS build the product of
        // transpose matrices, and read the result as if it were row-wise
        // again. In other words, we calculate (B^T A)^T, which is A^T B.

        const int m = src.n();
        const int n = this->n();
        const int k = this->m();
        const char *trans = "t";
        const char *notrans = "n";

        const number alpha = 1.;
        const number beta = (adding == true) ? 1. : 0.;

        // Use the BLAS function gemm for calculating the matrix-matrix
        // product.
        gemm(notrans, trans, &m, &n, &k, &alpha, &src(0,0), &m,
             &this->values[0], &n, &beta, &dst(0,0), &m);

        return;
      }

#endif

  const size_type m = n(), n = src.n(), l = this->m();

  // symmetric matrix if the two matrices are the same
  if (PointerComparison::equal(this, &src))
    for (size_type i=0; i<m; ++i)
      for (size_type j=i; j<m; ++j)
        {
          number2 add_value = 0.;
          for (size_type k=0; k<l; ++k)
            add_value += (number2)(*this)(k,i) * (number2)(*this)(k,j);
          if (adding)
            {
              dst(i,j) += add_value;
              if (i<j)
                dst(j,i) += add_value;
            }
          else
            dst(i,j) = dst(j,i) = add_value;
        }
  // arrange the loops in a way that we keep write operations low, (writing is
  // usually more costly than reading), even though we need to access the data
  // in src not in a contiguous way. However, we should usually end up in the
  // optimized gemm operation in case the matrix is big, so this shouldn't be
  // too bad.
  else
    for (size_type i=0; i<m; i++)
      for (size_type j=0; j<n; j++)
        {
          number2 add_value = adding ? dst(i,j) : 0.;
          for (size_type k=0; k<l; k++)
            add_value += (number2)(*this)(k,i) * (number2)(src(k,j));
          dst(i,j) = add_value;
        }
}



template <typename number>
template <typename number2>
void FullMatrix<number>::mTmult (FullMatrix<number2>       &dst,
                                 const FullMatrix<number2> &src,
                                 const bool                 adding) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));
  Assert (dst.n() == src.m(), ExcDimensionMismatch(dst.n(), src.m()));
  Assert (dst.m() == m(), ExcDimensionMismatch(m(), dst.m()));

  // see if we can use BLAS algorithms for this and if the type for 'number'
  // works for us (it is usually not efficient to use BLAS for very small
  // matrices):
#if defined(HAVE_DGEMM_) && defined (HAVE_SGEMM_)
  if ((types_are_equal<number,double>::value
       ||
       types_are_equal<number,float>::value)
      &&
      types_are_equal<number,number2>::value)
    if (this->n()*this->m()*src.m() > 300)
      {
        // In case we have the BLAS function gemm detected at configure, we
        // use that algorithm for matrix-matrix multiplication since it
        // provides better performance than the deal.II native function (it
        // uses cache and register blocking in order to access local data).
        //
        // Note that BLAS/LAPACK stores matrix elements column-wise (i.e., all
        // values in one column, then all in the next, etc.), whereas the
        // FullMatrix stores them row-wise.  We ignore that difference, and
        // give our row-wise data to BLAS, let BLAS build the product of
        // transpose matrices, and read the result as if it were row-wise
        // again. In other words, we calculate (B A^T)^T, which is AB^T.

        const int m = src.m();
        const int n = this->m();
        const int k = this->n();
        const char *notrans = "n";
        const char *trans = "t";

        const number alpha = 1.;
        const number beta = (adding == true) ? 1. : 0.;

        // Use the BLAS function gemm for calculating the matrix-matrix
        // product.
        gemm(trans, notrans, &m, &n, &k, &alpha, &src(0,0), &k,
             &this->values[0], &k, &beta, &dst(0,0), &m);

        return;
      }

#endif

  const size_type m = this->m(), n = src.m(), l = this->n();

  // symmetric matrix if the two matrices are the same
  if (PointerComparison::equal(this, &src))
    for (size_type i=0; i<m; ++i)
      for (size_type j=i; j<m; ++j)
        {
          number2 add_value = 0.;
          for (size_type k=0; k<l; ++k)
            add_value += (number2)(*this)(i,k) * (number2)(*this)(j,k);
          if (adding)
            {
              dst(i,j) += add_value;
              if (i<j)
                dst(j,i) += add_value;
            }
          else
            dst(i,j) = dst(j,i) = add_value;
        }
  else
    // arrange the loops in a way that we keep write operations low, (writing is
    // usually more costly than reading).
    for (size_type i=0; i<m; i++)
      for (size_type j=0; j<n; j++)
        {
          number2 add_value = adding ? dst(i,j) : 0.;
          for (size_type k=0; k<l; k++)
            add_value += (number2)(*this)(i,k) * (number2)(src(j,k));
          dst(i,j) = add_value;
        }
}



template <typename number>
template <typename number2>
void FullMatrix<number>::TmTmult (FullMatrix<number2>       &dst,
                                  const FullMatrix<number2> &src,
                                  const bool                 adding) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  Assert (m() == src.n(), ExcDimensionMismatch(m(), src.n()));
  Assert (n() == dst.m(), ExcDimensionMismatch(n(), dst.m()));
  Assert (src.m() == dst.n(), ExcDimensionMismatch(src.m(), dst.n()));


  // see if we can use BLAS algorithms for this and if the type for 'number'
  // works for us (it is usually not efficient to use BLAS for very small
  // matrices):
#if defined(HAVE_DGEMM_) && defined (HAVE_SGEMM_)
  if ((types_are_equal<number,double>::value
       ||
       types_are_equal<number,float>::value)
      &&
      types_are_equal<number,number2>::value)
    if (this->n()*this->m()*src.m() > 300)
      {
        // In case we have the BLAS function gemm detected at configure, we
        // use that algorithm for matrix-matrix multiplication since it
        // provides better performance than the deal.II native function (it
        // uses cache and register blocking in order to access local data).
        //
        // Note that BLAS/LAPACK stores matrix elements column-wise (i.e., all
        // values in one column, then all in the next, etc.), whereas the
        // FullMatrix stores them row-wise.  We ignore that difference, and
        // give our row-wise data to BLAS, let BLAS build the product of
        // transpose matrices, and read the result as if it were row-wise
        // again. In other words, we calculate (B A)^T, which is A^T B^T.

        const int m = src.n();
        const int n = this->n();
        const int k = this->m();
        const char *trans = "t";

        const number alpha = 1.;
        const number beta = (adding == true) ? 1. : 0.;

        // Use the BLAS function gemm for calculating the matrix-matrix
        // product.
        gemm(trans, trans, &m, &n, &k, &alpha, &src(0,0), &k,
             &this->values[0], &n, &beta, &dst(0,0), &m);

        return;
      }

#endif

  const size_type m = n(), n = src.m(), l = this->m();

  // arrange the loops in a way that we keep write operations low, (writing is
  // usually more costly than reading), even though we need to access the data
  // in the calling matrix in a non-contiguous way, possibly leading to cache
  // misses. However, we should usually end up in the optimized gemm operation
  // in case the matrix is big, so this shouldn't be too bad.
  for (size_type i=0; i<m; i++)
    for (size_type j=0; j<n; j++)
      {
        number2 add_value = adding ? dst(i,j) : 0.;
        for (size_type k=0; k<l; k++)
          add_value += (number2)(*this)(k,i) * (number2)(src(j,k));
        dst(i,j) = add_value;
      }
}


template <typename number>
void
FullMatrix<number>::triple_product(
  const FullMatrix<number> &A,
  const FullMatrix<number> &B,
  const FullMatrix<number> &D,
  const bool transpose_B,
  const bool transpose_D,
  const number scaling)
{
  if (transpose_B)
    {
      AssertDimension(B.m(), A.m());
      AssertDimension(B.n(), m());
    }
  else
    {
      AssertDimension(B.n(), A.m());
      AssertDimension(B.m(), m());
    }
  if (transpose_D)
    {
      AssertDimension(D.n(), A.n());
      AssertDimension(D.m(), n());
    }
  else
    {
      AssertDimension(D.m(), A.n());
      AssertDimension(D.n(), n());
    }

  // For all entries of the product
  // AD
  for (size_type i=0; i<A.m(); ++i)
    for (size_type j=0; j<n(); ++j)
      {
        // Compute the entry
        number ADij = 0.;
        if (transpose_D)
          for (size_type k=0; k<A.n(); ++k)
            ADij += A(i,k)*D(j,k);
        else
          for (size_type k=0; k<A.n(); ++k)
            ADij += A(i,k)*D(k,j);
        // And add it to this after
        // multiplying with the right
        // factor from B
        if (transpose_B)
          for (size_type k=0; k<m(); ++k)
            this->operator()(k,j) += scaling * ADij * B(i,k);
        else
          for (size_type k=0; k<m(); ++k)
            this->operator()(k,j) += scaling * ADij * B(k,i);
      }
}


template <typename number>
template <typename number2>
number2
FullMatrix<number>::matrix_norm_square (const Vector<number2> &v) const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert(m() == v.size(), ExcDimensionMismatch(m(),v.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(),v.size()));

  number2 sum = 0.;
  const size_type n_rows = m();
  const number *val_ptr = &this->values[0];
  const number2 *v_ptr;

  for (size_type row=0; row<n_rows; ++row)
    {
      number s = 0.;
      const number *const val_end_of_row = val_ptr+n_rows;
      v_ptr = v.begin();
      while (val_ptr != val_end_of_row)
        s += number(*val_ptr++) * number(*v_ptr++);

      sum += s * number(numbers::NumberTraits<number2>::conjugate(v(row)));
    }

  return sum;
}


template <typename number>
template <typename number2>
number2
FullMatrix<number>::matrix_scalar_product (const Vector<number2> &u,
                                           const Vector<number2> &v) const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert(m() == u.size(), ExcDimensionMismatch(m(),v.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(),v.size()));

  number2 sum = 0.;
  const size_type n_rows = m();
  const size_type n_cols = n();
  const number *val_ptr = &this->values[0];
  const number2 *v_ptr;

  for (size_type row=0; row<n_rows; ++row)
    {
      number s = 0.;
      const number *const val_end_of_row = val_ptr+n_cols;
      v_ptr = v.begin();
      while (val_ptr != val_end_of_row)
        s += number(*val_ptr++) * number(*v_ptr++);

      sum += s * number(u(row));
    }

  return sum;
}



template <typename number>
void
FullMatrix<number>::symmetrize ()
{
  Assert (m() == n(), ExcNotQuadratic());

  const size_type N = m();
  for (size_type i=0; i<N; ++i)
    for (size_type j=i+1; j<N; ++j)
      {
        const number t = ((*this)(i,j) + (*this)(j,i)) / number(2.);
        (*this)(i,j) = (*this)(j,i) = t;
      };
}


template <typename number>
typename FullMatrix<number>::real_type
FullMatrix<number>::l1_norm () const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  real_type sum=0, max=0;
  const size_type n_rows = m(), n_cols = n();

  for (size_type col=0; col<n_cols; ++col)
    {
      sum=0;
      for (size_type row=0; row<n_rows; ++row)
        sum += std::abs((*this)(row,col));
      if (sum > max)
        max = sum;
    }
  return max;
}



template <typename number>
typename FullMatrix<number>::real_type
FullMatrix<number>::linfty_norm () const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  real_type sum=0, max=0;
  const size_type n_rows = m(), n_cols = n();

  for (size_type row=0; row<n_rows; ++row)
    {
      sum=0;
      for (size_type col=0; col<n_cols; ++col)
        sum += std::abs((*this)(row,col));
      if (sum > max)
        max = sum;
    }
  return max;
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::add (const number               a,
                         const FullMatrix<number2> &A)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (m() == A.m(), ExcDimensionMismatch(m(), A.m()));
  Assert (n() == A.n(), ExcDimensionMismatch(n(), A.n()));

  for (size_type i=0; i<m(); ++i)
    for (size_type j=0; j<n(); ++j)
      (*this)(i,j) += a * number(A(i,j));
}


template <typename number>
template <typename number2>
void
FullMatrix<number>::add (const number               a,
                         const FullMatrix<number2> &A,
                         const number               b,
                         const FullMatrix<number2> &B)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (m() == A.m(), ExcDimensionMismatch(m(), A.m()));
  Assert (n() == A.n(), ExcDimensionMismatch(n(), A.n()));
  Assert (m() == B.m(), ExcDimensionMismatch(m(), B.m()));
  Assert (n() == B.n(), ExcDimensionMismatch(n(), B.n()));

  for (size_type i=0; i<m(); ++i)
    for (size_type j=0; j<n(); ++j)
      (*this)(i,j) += a * number(A(i,j)) + b * number(B(i,j));
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::add (const number               a,
                         const FullMatrix<number2> &A,
                         const number               b,
                         const FullMatrix<number2> &B,
                         const number               c,
                         const FullMatrix<number2> &C)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (m() == A.m(), ExcDimensionMismatch(m(), A.m()));
  Assert (n() == A.n(), ExcDimensionMismatch(n(), A.n()));
  Assert (m() == B.m(), ExcDimensionMismatch(m(), B.m()));
  Assert (n() == B.n(), ExcDimensionMismatch(n(), B.n()));
  Assert (m() == C.m(), ExcDimensionMismatch(m(), C.m()));
  Assert (n() == C.n(), ExcDimensionMismatch(n(), C.n()));


  for (size_type i=0; i<m(); ++i)
    for (size_type j=0; j<n(); ++j)
      (*this)(i,j) += a * number(A(i,j)) +
                      b * number(B(i,j)) +
                      c * number(C(i,j));
}



template <typename number>
template <typename number2>
void FullMatrix<number>::add (const FullMatrix<number2> &src,
                              const number factor,
                              const size_type dst_offset_i,
                              const size_type dst_offset_j,
                              const size_type src_offset_i,
                              const size_type src_offset_j)
{
  Assert (dst_offset_i < m(),
          ExcIndexRange (dst_offset_i, 0, m()));
  Assert (dst_offset_j < n(),
          ExcIndexRange (dst_offset_j, 0, n()));
  Assert (src_offset_i < src.m(),
          ExcIndexRange (src_offset_i, 0, src.m()));
  Assert (src_offset_j < src.n(),
          ExcIndexRange (src_offset_j, 0, src.n()));

  // Compute maximal size of copied block
  const size_type rows = std::min (m() - dst_offset_i, src.m() - src_offset_i);
  const size_type cols = std::min (n() - dst_offset_j, src.n() - src_offset_j);

  for (size_type i=0; i<rows ; ++i)
    for (size_type j=0; j<cols ; ++j)
      (*this)(dst_offset_i+i,dst_offset_j+j)
      += factor * number(src(src_offset_i+i,src_offset_j+j));
}



template <typename number>
template <typename number2>
void FullMatrix<number>::Tadd (const FullMatrix<number2> &src,
                               const number factor,
                               const size_type dst_offset_i,
                               const size_type dst_offset_j,
                               const size_type src_offset_i,
                               const size_type src_offset_j)
{
  Assert (dst_offset_i < m(),
          ExcIndexRange (dst_offset_i, 0, m()));
  Assert (dst_offset_j < n(),
          ExcIndexRange (dst_offset_j, 0, n()));
  Assert (src_offset_i < src.n(),
          ExcIndexRange (src_offset_i, 0, src.n()));
  Assert (src_offset_j < src.m(),
          ExcIndexRange (src_offset_j, 0, src.m()));

  // Compute maximal size of copied block
  const size_type rows = std::min (m() - dst_offset_i, src.n() - src_offset_j);
  const size_type cols = std::min (n() - dst_offset_j,
                                   src.m() - src_offset_i);


  for (size_type i=0; i<rows ; ++i)
    for (size_type j=0; j<cols ; ++j)
      (*this)(dst_offset_i+i,dst_offset_j+j)
      += factor * number(src(src_offset_i+j,src_offset_j+i));
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::Tadd (const number a,
                          const FullMatrix<number2> &A)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (m() == n(),     ExcNotQuadratic());
  Assert (m() == A.m(), ExcDimensionMismatch(m(), A.m()));
  Assert (n() == A.n(), ExcDimensionMismatch(n(), A.n()));

  for (size_type i=0; i<n(); ++i)
    for (size_type j=0; j<m(); ++j)
      (*this)(i,j) += a * number(A(j,i));
}


template <typename number>
bool
FullMatrix<number>::operator == (const FullMatrix<number> &M) const
{
  // simply pass down to the base class
  return Table<2,number>::operator==(M);
}


template <typename number>
number
FullMatrix<number>::determinant () const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (this->n_cols() == this->n_rows(),
          ExcDimensionMismatch(this->n_cols(), this->n_rows()));

  switch (this->n_cols())
    {
    case 1:
      return (*this)(0,0);
    case 2:
      return (*this)(0,0)*(*this)(1,1) - (*this)(1,0)*(*this)(0,1);
    case 3:
      return  ((*this)(0,0)*(*this)(1,1)*(*this)(2,2)
               -(*this)(0,0)*(*this)(1,2)*(*this)(2,1)
               -(*this)(1,0)*(*this)(0,1)*(*this)(2,2)
               +(*this)(1,0)*(*this)(0,2)*(*this)(2,1)
               +(*this)(2,0)*(*this)(0,1)*(*this)(1,2)
               -(*this)(2,0)*(*this)(0,2)*(*this)(1,1));
    default:
      Assert (false, ExcNotImplemented());
      return 0;
    };
}



template <typename number>
number
FullMatrix<number>::trace () const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (this->n_cols() == this->n_rows(),
          ExcDimensionMismatch(this->n_cols(), this->n_rows()));

  number tr = 0;
  for (size_type i=0; i<this->n_rows(); ++i)
    tr += (*this)(i,i);

  return tr;
}



template <typename number>
typename FullMatrix<number>::real_type
FullMatrix<number>::frobenius_norm () const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  real_type s = 0.;
  for (size_type i=0; i<this->n_rows()*this->n_cols(); ++i)
    s += numbers::NumberTraits<number>::abs_square(this->values[i]);
  return std::sqrt(s);
}



template <typename number>
typename FullMatrix<number>::real_type
FullMatrix<number>::relative_symmetry_norm2 () const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  real_type s = 0.;
  real_type a = 0.;
  for (size_type i=0; i<this->n_rows(); ++i)
    for (size_type j=0; j<this->n_cols(); ++j)
      {
        const number x_ij = (*this)(i,j);
        const number x_ji = (*this)(j,i);

        a += numbers::NumberTraits<number>::abs_square(x_ij-x_ji);
        s += numbers::NumberTraits<number>::abs_square(x_ij);
      }

  if (s!=0.)
    return std::sqrt(a)/std::sqrt(s);
  return 0;
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::invert (const FullMatrix<number2> &M)
{
  Assert (!this->empty(), ExcEmptyMatrix());

  Assert (this->n_cols() == this->n_rows(),
          ExcNotQuadratic());
  Assert (this->n_cols() == M.n_cols(),
          ExcDimensionMismatch(this->n_cols(), M.n_cols()));
  Assert (this->n_rows() == M.n_rows(),
          ExcDimensionMismatch(this->n_rows(), M.n_rows()));

  if (PointerComparison::equal(&M, this))
    {
      // avoid overwriting source
      // by destination matrix:
      FullMatrix<number2> M2 = M;
      invert(M2);
    }
  else
    switch (this->n_cols())
      {
      case 1:
        (*this)(0,0) = number2(1.0)/M(0,0);
        return;
      case 2:
        // this is Maple output,
        // thus a bit unstructured
      {
        const number2 t4 = number2(1.0)/(M(0,0)*M(1,1)-M(0,1)*M(1,0));
        (*this)(0,0) = M(1,1)*t4;
        (*this)(0,1) = -M(0,1)*t4;
        (*this)(1,0) = -M(1,0)*t4;
        (*this)(1,1) = M(0,0)*t4;
        return;
      };

      case 3:
      {
        const number2 t4 = M(0,0)*M(1,1),
                      t6 = M(0,0)*M(1,2),
                      t8 = M(0,1)*M(1,0),
                      t00 = M(0,2)*M(1,0),
                      t01 = M(0,1)*M(2,0),
                      t04 = M(0,2)*M(2,0),
                      t07 = number2(1.0)/(t4*M(2,2)-t6*M(2,1)-t8*M(2,2)+
                                          t00*M(2,1)+t01*M(1,2)-t04*M(1,1));
        (*this)(0,0) = (M(1,1)*M(2,2)-M(1,2)*M(2,1))*t07;
        (*this)(0,1) = -(M(0,1)*M(2,2)-M(0,2)*M(2,1))*t07;
        (*this)(0,2) = -(-M(0,1)*M(1,2)+M(0,2)*M(1,1))*t07;
        (*this)(1,0) = -(M(1,0)*M(2,2)-M(1,2)*M(2,0))*t07;
        (*this)(1,1) = (M(0,0)*M(2,2)-t04)*t07;
        (*this)(1,2) = -(t6-t00)*t07;
        (*this)(2,0) = -(-M(1,0)*M(2,1)+M(1,1)*M(2,0))*t07;
        (*this)(2,1) = -(M(0,0)*M(2,1)-t01)*t07;
        (*this)(2,2) = (t4-t8)*t07;
        return;
      };

      case 4:
      {
        // with (linalg);
        // a:=matrix(4,4);
        // evalm(a);
        // ai:=inverse(a);
        // readlib(C);
        // C(ai,optimized,filename=x4);

        const number2 t14 = M(0,0)*M(1,1);
        const number2 t15 = M(2,2)*M(3,3);
        const number2 t17 = M(2,3)*M(3,2);
        const number2 t19 = M(0,0)*M(2,1);
        const number2 t20 = M(1,2)*M(3,3);
        const number2 t22 = M(1,3)*M(3,2);
        const number2 t24 = M(0,0)*M(3,1);
        const number2 t25 = M(1,2)*M(2,3);
        const number2 t27 = M(1,3)*M(2,2);
        const number2 t29 = M(1,0)*M(0,1);
        const number2 t32 = M(1,0)*M(2,1);
        const number2 t33 = M(0,2)*M(3,3);
        const number2 t35 = M(0,3)*M(3,2);
        const number2 t37 = M(1,0)*M(3,1);
        const number2 t38 = M(0,2)*M(2,3);
        const number2 t40 = M(0,3)*M(2,2);
        const number2 t42 = t14*t15-t14*t17-t19*t20+t19*t22+
                            t24*t25-t24*t27-t29*t15+t29*t17+
                            t32*t33-t32*t35-t37*t38+t37*t40;
        const number2 t43 = M(2,0)*M(0,1);
        const number2 t46 = M(2,0)*M(1,1);
        const number2 t49 = M(2,0)*M(3,1);
        const number2 t50 = M(0,2)*M(1,3);
        const number2 t52 = M(0,3)*M(1,2);
        const number2 t54 = M(3,0)*M(0,1);
        const number2 t57 = M(3,0)*M(1,1);
        const number2 t60 = M(3,0)*M(2,1);
        const number2 t63 = t43*t20-t43*t22-t46*t33+t46*t35+
                            t49*t50-t49*t52-t54*t25+t54*t27+
                            t57*t38-t57*t40-t60*t50+t60*t52;
        const number2 t65 = number2(1.)/(t42+t63);
        const number2 t71 = M(0,2)*M(2,1);
        const number2 t73 = M(0,3)*M(2,1);
        const number2 t75 = M(0,2)*M(3,1);
        const number2 t77 = M(0,3)*M(3,1);
        const number2 t81 = M(0,1)*M(1,2);
        const number2 t83 = M(0,1)*M(1,3);
        const number2 t85 = M(0,2)*M(1,1);
        const number2 t87 = M(0,3)*M(1,1);
        const number2 t101 = M(1,0)*M(2,2);
        const number2 t103 = M(1,0)*M(2,3);
        const number2 t105 = M(2,0)*M(1,2);
        const number2 t107 = M(2,0)*M(1,3);
        const number2 t109 = M(3,0)*M(1,2);
        const number2 t111 = M(3,0)*M(1,3);
        const number2 t115 = M(0,0)*M(2,2);
        const number2 t117 = M(0,0)*M(2,3);
        const number2 t119 = M(2,0)*M(0,2);
        const number2 t121 = M(2,0)*M(0,3);
        const number2 t123 = M(3,0)*M(0,2);
        const number2 t125 = M(3,0)*M(0,3);
        const number2 t129 = M(0,0)*M(1,2);
        const number2 t131 = M(0,0)*M(1,3);
        const number2 t133 = M(1,0)*M(0,2);
        const number2 t135 = M(1,0)*M(0,3);
        (*this)(0,0) = (M(1,1)*M(2,2)*M(3,3)-M(1,1)*M(2,3)*M(3,2)-
                        M(2,1)*M(1,2)*M(3,3)+M(2,1)*M(1,3)*M(3,2)+
                        M(3,1)*M(1,2)*M(2,3)-M(3,1)*M(1,3)*M(2,2))*t65;
        (*this)(0,1) = -(M(0,1)*M(2,2)*M(3,3)-M(0,1)*M(2,3)*M(3,2)-
                         t71*M(3,3)+t73*M(3,2)+t75*M(2,3)-t77*M(2,2))*t65;
        (*this)(0,2) = (t81*M(3,3)-t83*M(3,2)-t85*M(3,3)+t87*M(3,2)+
                        t75*M(1,3)-t77*M(1,2))*t65;
        (*this)(0,3) = -(t81*M(2,3)-t83*M(2,2)-t85*M(2,3)+t87*M(2,2)+
                         t71*M(1,3)-t73*M(1,2))*t65;
        (*this)(1,0) = -(t101*M(3,3)-t103*M(3,2)-t105*M(3,3)+t107*M(3,2)+
                         t109*M(2,3)-t111*M(2,2))*t65;
        (*this)(1,1) = (t115*M(3,3)-t117*M(3,2)-t119*M(3,3)+t121*M(3,2)+
                        t123*M(2,3)-t125*M(2,2))*t65;
        (*this)(1,2) = -(t129*M(3,3)-t131*M(3,2)-t133*M(3,3)+t135*M(3,2)+
                         t123*M(1,3)-t125*M(1,2))*t65;
        (*this)(1,3) = (t129*M(2,3)-t131*M(2,2)-t133*M(2,3)+t135*M(2,2)+
                        t119*M(1,3)-t121*M(1,2))*t65;
        (*this)(2,0) = (t32*M(3,3)-t103*M(3,1)-t46*M(3,3)+t107*M(3,1)+
                        t57*M(2,3)-t111*M(2,1))*t65;
        (*this)(2,1) = -(t19*M(3,3)-t117*M(3,1)-t43*M(3,3)+t121*M(3,1)+
                         t54*M(2,3)-t125*M(2,1))*t65;
        (*this)(2,2) = (t14*M(3,3)-t131*M(3,1)-t29*M(3,3)+t135*M(3,1)+
                        t54*M(1,3)-t125*M(1,1))*t65;
        (*this)(2,3) = -(t14*M(2,3)-t131*M(2,1)-t29*M(2,3)+t135*M(2,1)+
                         t43*M(1,3)-t121*M(1,1))*t65;
        (*this)(3,0) = -(t32*M(3,2)-t101*M(3,1)-t46*M(3,2)+t105*M(3,1)+
                         t57*M(2,2)-t109*M(2,1))*t65;
        (*this)(3,1) = (t19*M(3,2)-t115*M(3,1)-t43*M(3,2)+t119*M(3,1)+
                        t54*M(2,2)-t123*M(2,1))*t65;
        (*this)(3,2) = -(t14*M(3,2)-t129*M(3,1)-t29*M(3,2)+t133*M(3,1)+
                         t54*M(1,2)-t123*M(1,1))*t65;
        (*this)(3,3) = (t14*M(2,2)-t129*M(2,1)-t29*M(2,2)+t133*M(2,1)+
                        t43*M(1,2)-t119*M(1,1))*t65;

        break;
      }


      default:
        // if no inversion is
        // hardcoded, fall back
        // to use the
        // Gauss-Jordan algorithm
        *this = M;
        gauss_jordan();
      };
}


template <typename number>
template <typename number2>
void
FullMatrix<number>::cholesky (const FullMatrix<number2> &A)
{
  Assert (!A.empty(), ExcEmptyMatrix());
  Assert (A.n() == A.m(),
          ExcNotQuadratic());
  // Matrix must be symmetric.
  Assert(A.relative_symmetry_norm2() < 1.0e-10, ExcMessage("A must be symmetric."));

  if (PointerComparison::equal(&A, this))
    {
      // avoid overwriting source
      // by destination matrix:
      FullMatrix<number2> A2 = A;
      cholesky(A2);
    }
  else
    {
      /* reinit *this to 0 */
      this->reinit(A.m(), A.n());

      double SLik2 = 0.0, SLikLjk = 0.0;
      for (size_type i=0; i< this->n_cols(); i++)
        {
          SLik2 = 0.0;
          for (size_type j = 0; j < i; j++)
            {
              SLikLjk = 0.0;
              for (size_type k =0; k<j; k++)
                {
                  SLikLjk += (*this)(i,k)*(*this)(j,k);
                };
              (*this)(i,j) = (1./(*this)(j,j))*(A(i,j) - SLikLjk);
              SLik2 += (*this)(i,j)*(*this)(i,j);
            }
          AssertThrow (A(i,i) - SLik2 >= 0,
                       ExcMatrixNotPositiveDefinite());

          (*this)(i,i) = std::sqrt(A(i,i) - SLik2);
        }
    }
}


template <typename number>
template <typename number2>
void
FullMatrix<number>::outer_product (const Vector<number2> &V,
                                   const Vector<number2> &W)
{
  Assert (V.size() == W.size(), ExcMessage("Vectors V, W must be the same size."));
  this->reinit(V.size(), V.size());

  for (size_type i = 0; i<this->n(); i++)
    {
      for (size_type j = 0; j< this->n(); j++)
        {
          (*this)(i,j) = V(i)*W(j);
        }
    }
}


template <typename number>
template <typename number2>
void
FullMatrix<number>::left_invert (const FullMatrix<number2> &A)
{
  Assert (!A.empty(), ExcEmptyMatrix());
  Assert(A.m()>A.n(), ExcDimensionMismatch(A.m(), A.n()));
  Assert(this->m()==A.n(), ExcDimensionMismatch(this->m(), A.n()));
  Assert(this->n()==A.m(), ExcDimensionMismatch(this->n(), A.m()));

  FullMatrix<number2> A_t(A.n(),A.m());
  FullMatrix<number2> A_t_times_A(A.n(),A.n());
  FullMatrix<number2> A_t_times_A_inv(A.n(),A.n());
  FullMatrix<number2> left_inv(A.n(),A.m());

  A_t.Tadd(A,1);
  A_t.mmult(A_t_times_A,A);
  if (number(A_t_times_A.determinant())==number(0))
    Assert(false, ExcSingular())
    else
      {
        A_t_times_A_inv.invert(A_t_times_A);
        A_t_times_A_inv.mmult(left_inv,A_t);

        *this=left_inv;
      }
}

template <typename number>
template <typename number2>
void
FullMatrix<number>::right_invert (const FullMatrix<number2> &A)
{
  Assert (!A.empty(), ExcEmptyMatrix());
  Assert(A.n()>A.m(), ExcDimensionMismatch(A.n(), A.m()));
  Assert(this->m()==A.n(), ExcDimensionMismatch(this->m(), A.n()));
  Assert(this->n()==A.m(), ExcDimensionMismatch(this->n(), A.m()));

  FullMatrix<number> A_t(A.n(),A.m());
  FullMatrix<number> A_times_A_t(A.m(),A.m());
  FullMatrix<number> A_times_A_t_inv(A.m(),A.m());
  FullMatrix<number> right_inv(A.n(),A.m());

  A_t.Tadd(A,1);
  A.mmult(A_times_A_t,A_t);
  if (number(A_times_A_t.determinant())==number(0))
    Assert(false, ExcSingular())
    else
      {
        A_times_A_t_inv.invert(A_times_A_t);
        A_t.mmult(right_inv,A_times_A_t_inv);

        *this=right_inv;
      }
}


template <typename number>
template <int dim>
void
FullMatrix<number>::copy_from (const Tensor<2,dim> &T,
                               const size_type src_r_i,
                               const size_type src_r_j,
                               const size_type src_c_i,
                               const size_type src_c_j,
                               const size_type dst_r,
                               const size_type dst_c)
{

  Assert (!this->empty(), ExcEmptyMatrix());
  Assert(this->m()-dst_r>src_r_j-src_r_i,
         ExcIndexRange(this->m()-dst_r,0,src_r_j-src_r_i));
  Assert(this->n()-dst_c>src_c_j-src_c_i,
         ExcIndexRange(this->n()-dst_c,0,src_c_j-src_c_i));
  Assert(dim>src_r_j, ExcIndexRange(dim,0,src_r_j));
  Assert(dim>src_c_j, ExcIndexRange(dim,0,src_r_j));
  Assert(src_r_j>=src_r_i, ExcIndexRange(src_r_j,0,src_r_i));
  Assert(src_c_j>=src_c_i, ExcIndexRange(src_r_j,0,src_r_i));

  for (size_type i=0; i<src_r_j-src_r_i+1; i++)
    for (size_type j=0; j<src_c_j-src_c_i+1; j++)
      (*this)(i+dst_r,j+dst_c) = number(T[i+src_r_i][j+src_c_i]);

}


template <typename number>
template <int dim>
void
FullMatrix<number>::copy_to (Tensor<2,dim> &T,
                             const size_type src_r_i,
                             const size_type src_r_j,
                             const size_type src_c_i,
                             const size_type src_c_j,
                             const size_type dst_r,
                             const size_type dst_c) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  Assert(dim-dst_r>src_r_j-src_r_i,
         ExcIndexRange(dim-dst_r,0,src_r_j-src_r_i));
  Assert(dim-dst_c>src_c_j-src_c_i,
         ExcIndexRange(dim-dst_c,0,src_c_j-src_c_i));
  Assert(this->m()>src_r_j, ExcIndexRange(dim,0,src_r_j));
  Assert(this->n()>src_c_j, ExcIndexRange(dim,0,src_r_j));
  Assert(src_r_j>=src_r_i, ExcIndexRange(src_r_j,0,src_r_i));
  Assert(src_c_j>=src_c_i, ExcIndexRange(src_r_j,0,src_r_i));


  for (size_type i=0; i<src_r_j-src_r_i+1; i++)
    for (size_type j=0; j<src_c_j-src_c_i+1; j++)
      T[i+dst_r][j+dst_c] = double ((*this)(i+src_r_i,j+src_c_i));
}



template <typename number>
template <typename somenumber>
void
FullMatrix<number>::precondition_Jacobi (Vector<somenumber>       &dst,
                                         const Vector<somenumber> &src,
                                         const number              om) const
{
  Assert (m() == n(), ExcNotQuadratic());
  Assert (dst.size() == n(), ExcDimensionMismatch (dst.size(), n()));
  Assert (src.size() == n(), ExcDimensionMismatch (src.size(), n()));

  const size_t n = src.size();
  somenumber       *dst_ptr = dst.begin();
  const somenumber *src_ptr = src.begin();

  for (size_type i=0; i<n; ++i, ++dst_ptr, ++src_ptr)
    *dst_ptr = somenumber(om) **src_ptr / somenumber((*this)(i,i));
}



template <typename number>
void
FullMatrix<number>::print_formatted (
  std::ostream       &out,
  const unsigned int  precision,
  const bool          scientific,
  const unsigned int  width_,
  const char         *zero_string,
  const double        denominator,
  const double        threshold) const
{
  unsigned int width = width_;

  Assert ((!this->empty()) || (this->n_cols()+this->n_rows()==0),
          ExcInternalError());

  // set output format, but store old
  // state
  std::ios::fmtflags old_flags = out.flags();
  unsigned int old_precision = out.precision (precision);

  if (scientific)
    {
      out.setf (std::ios::scientific, std::ios::floatfield);
      if (!width)
        width = precision+7;
    }
  else
    {
      out.setf (std::ios::fixed, std::ios::floatfield);
      if (!width)
        width = precision+2;
    }

  for (size_type i=0; i<m(); ++i)
    {
      for (size_type j=0; j<n(); ++j)
        if (std::abs((*this)(i,j)) > threshold)
          out << std::setw(width)
              << (*this)(i,j) * number(denominator) << ' ';
        else
          out << std::setw(width) << zero_string << ' ';
      out << std::endl;
    };

  AssertThrow (out, ExcIO());
  // reset output format
  out.flags (old_flags);
  out.precision(old_precision);
}


template <typename number>
void
FullMatrix<number>::gauss_jordan ()
{
  Assert (!this->empty(), ExcEmptyMatrix());
  Assert (this->n_cols() == this->n_rows(), ExcNotQuadratic());

  // see if we can use Lapack algorithms
  // for this and if the type for 'number'
  // works for us (it is usually not
  // efficient to use Lapack for very small
  // matrices):
#if defined(HAVE_DGETRF_) && defined (HAVE_SGETRF_) && \
    defined(HAVE_DGETRI_) && defined (HAVE_SGETRI_)
  if (types_are_equal<number,double>::value
      ||
      types_are_equal<number,float>::value)
    if (this->n_cols() > 15)
      {
        // In case we have the LAPACK functions
        // getrf and getri detected at configure,
        // we use these algorithms for inversion
        // since they provide better performance
        // than the deal.II native functions.
        //
        // Note that BLAS/LAPACK stores matrix
        // elements column-wise (i.e., all values in
        // one column, then all in the next, etc.),
        // whereas the FullMatrix stores them
        // row-wise.
        // We ignore that difference, and give our
        // row-wise data to LAPACK,
        // let LAPACK build the inverse of the
        // transpose matrix, and read the result as
        // if it were row-wise again. In other words,
        // we just got ((A^T)^{-1})^T, which is
        // A^{-1}.

        const int nn = this->n();

        // workspace for permutations
        std::vector<int> ipiv(nn);
        int info;

        // Use the LAPACK function getrf for
        // calculating the LU factorization.
        getrf(&nn, &nn, &this->values[0], &nn, &ipiv[0], &info);

        Assert(info >= 0, ExcInternalError());
        Assert(info == 0, LACExceptions::ExcSingular());

        // scratch array
        std::vector<number> inv_work (nn);

        // Use the LAPACK function getri for
        // calculating the actual inverse using
        // the LU factorization.
        getri(&nn, &this->values[0], &nn, &ipiv[0], &inv_work[0], &nn, &info);

        Assert(info >= 0, ExcInternalError());
        Assert(info == 0, LACExceptions::ExcSingular());

        return;
      }

#endif

  // otherwise do it by hand. use the
  // Gauss-Jordan-Algorithmus from
  // Stoer & Bulirsch I (4th Edition)
  // p. 153
  const size_type N = n();

  // first get an estimate of the
  // size of the elements of this
  // matrix, for later checks whether
  // the pivot element is large
  // enough, or whether we have to
  // fear that the matrix is not
  // regular
  double diagonal_sum = 0;
  for (size_type i=0; i<N; ++i)
    diagonal_sum += std::abs((*this)(i,i));
  const double typical_diagonal_element = diagonal_sum/N;

  // initialize the array that holds
  // the permutations that we find
  // during pivot search
  std::vector<size_type> p(N);
  for (size_type i=0; i<N; ++i)
    p[i] = i;

  for (size_type j=0; j<N; ++j)
    {
      // pivot search: search that
      // part of the line on and
      // right of the diagonal for
      // the largest element
      real_type max = std::abs((*this)(j,j));
      size_type r   = j;
      for (size_type i=j+1; i<N; ++i)
        {
          if (std::abs((*this)(i,j)) > max)
            {
              max = std::abs((*this)(i,j));
              r = i;
            }
        }
      // check whether the pivot is
      // too small
      Assert(max > 1.e-16*typical_diagonal_element,
             ExcNotRegular(max));

      // row interchange
      if (r>j)
        {
          for (size_type k=0; k<N; ++k)
            std::swap ((*this)(j,k), (*this)(r,k));

          std::swap (p[j], p[r]);
        }

      // transformation
      const number hr = number(1.)/(*this)(j,j);
      (*this)(j,j) = hr;
      for (size_type k=0; k<N; ++k)
        {
          if (k==j) continue;
          for (size_type i=0; i<N; ++i)
            {
              if (i==j) continue;
              (*this)(i,k) -= (*this)(i,j)*(*this)(j,k)*hr;
            }
        }
      for (size_type i=0; i<N; ++i)
        {
          (*this)(i,j) *= hr;
          (*this)(j,i) *= -hr;
        }
      (*this)(j,j) = hr;
    }
  // column interchange
  std::vector<number> hv(N);
  for (size_type i=0; i<N; ++i)
    {
      for (size_type k=0; k<N; ++k)
        hv[p[k]] = (*this)(i,k);
      for (size_type k=0; k<N; ++k)
        (*this)(i,k) = hv[k];
    }
}



template <typename number>
std::size_t
FullMatrix<number>::memory_consumption () const
{
  return (sizeof(*this) - sizeof (Table<2,number>)
          +
          Table<2,number>::memory_consumption());
}


DEAL_II_NAMESPACE_CLOSE

#endif
