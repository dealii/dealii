//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2008, 2010 by the deal.II authors
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
#include <lac/block_vector.h>

#include <iostream>
#include <iomanip>

DEAL_II_NAMESPACE_OPEN

using namespace LAPACKSupport;

template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix(const unsigned int n)
		:
		TransposeTable<number> (n,n),
		state(matrix)
{}


template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix(
  const unsigned int m,
  const unsigned int n)
		:
		TransposeTable<number> (m,n),
		state(matrix)
{}


template <typename number>
LAPACKFullMatrix<number>::LAPACKFullMatrix(const LAPACKFullMatrix &M)
		:
		TransposeTable<number> (M),
		state(matrix)
{}


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
    std::fill (this->values.begin(), this->values.end(), number());

  state = LAPACKSupport::matrix;
  return *this;
}


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
  const number null = 0.;
  
  switch (state)
    {
      case matrix:
      case inverse_matrix:
      {	
	AssertDimension(v.size(), this->n_cols());
	AssertDimension(w.size(), this->n_rows());
	
	gemv("N", &mm, &nn, &alpha, this->data(), &mm, v.val, &one, &beta, w.val, &one);
	break;
      }
      case svd:
      {
	AssertDimension(v.size(), this->n_cols());
	AssertDimension(w.size(), this->n_rows());
					 // Compute V^T v
	work.resize(std::max(mm,nn));
	gemv("N", &nn, &nn, &alpha, svd_vt->data(), &nn, v.val, &one, &null, &work[0], &one);
					 // Multiply by singular values
	for (unsigned int i=0;i<wr.size();++i)
	  work[i] *= wr[i];
					 // Multiply with U
	gemv("N", &mm, &mm, &alpha, svd_u->data(), &mm, &work[0], &one, &beta, w.val, &one);
	break;
      }
      case inverse_svd:
      {
	AssertDimension(w.size(), this->n_cols());
	AssertDimension(v.size(), this->n_rows());
					 // Compute U^T v
	work.resize(std::max(mm,nn));
	gemv("T", &mm, &mm, &alpha, svd_u->data(), &mm, v.val, &one, &null, &work[0], &one);
					 // Multiply by singular values
	for (unsigned int i=0;i<wr.size();++i)
	  work[i] *= wr[i];
					 // Multiply with V
	gemv("T", &nn, &nn, &alpha, svd_vt->data(), &nn, &work[0], &one, &beta, w.val, &one);
	break;
      }
      default:
	    Assert (false, ExcState(state));
    }
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
  const number null = 0.;
  
  switch (state)
    {
      case matrix:
      case inverse_matrix:
      {	
	AssertDimension(w.size(), this->n_cols());
	AssertDimension(v.size(), this->n_rows());
	
	gemv("T", &mm, &nn, &alpha, this->data(), &mm, v.val, &one, &beta, w.val, &one);
	break;
      }
      case svd:
      {
	AssertDimension(w.size(), this->n_cols());
	AssertDimension(v.size(), this->n_rows());
	
					 // Compute U^T v
	work.resize(std::max(mm,nn));
	gemv("T", &mm, &mm, &alpha, svd_u->data(), &mm, v.val, &one, &null, &work[0], &one);
					 // Multiply by singular values
	for (unsigned int i=0;i<wr.size();++i)
	  work[i] *= wr[i];
					 // Multiply with V
	gemv("T", &nn, &nn, &alpha, svd_vt->data(), &nn, &work[0], &one, &beta, w.val, &one);
	break;
      case inverse_svd:
      {
	AssertDimension(v.size(), this->n_cols());
	AssertDimension(w.size(), this->n_rows());
  
					 // Compute V^T v
	work.resize(std::max(mm,nn));
	gemv("N", &nn, &nn, &alpha, svd_vt->data(), &nn, v.val, &one, &null, &work[0], &one);
					 // Multiply by singular values
	for (unsigned int i=0;i<wr.size();++i)
	  work[i] *= wr[i];
					 // Multiply with U
	gemv("N", &mm, &mm, &alpha, svd_u->data(), &mm, &work[0], &one, &beta, w.val, &one);
	break;
      }
      }
      default:
	    Assert (false, ExcState(state));
    }
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_lu_factorization()
{
  Assert(state == matrix, ExcState(state));
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  number* values = const_cast<number*> (this->data());
  ipiv.resize(mm);
  int info;
  getrf(&mm, &nn, values, &mm, &ipiv[0], &info);

  Assert(info >= 0, ExcInternalError());
  Assert(info == 0, LACExceptions::ExcSingular());

  state = lu;
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_svd()
{
  Assert(state == matrix, ExcState(state));
  state = LAPACKSupport::unusable;
  
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  number* values = const_cast<number*> (this->data());
  wr.resize(std::max(mm,nn));
  std::fill(wr.begin(), wr.end(), 0.);
  ipiv.resize(8*mm);
  
  svd_u.reset (new LAPACKFullMatrix<number>(mm,mm));
  svd_vt.reset (new LAPACKFullMatrix<number>(nn,nn));
  number* mu  = const_cast<number*> (svd_u->data());
  number* mvt = const_cast<number*> (svd_vt->data());
  int info;

				   // see comment on this #if
				   // below. Another reason to love Petsc
#ifndef DEAL_II_LIBLAPACK_NOQUERYMODE

				   // First determine optimal
				   // workspace size
  work.resize(1);
  int lwork = -1;
  gesdd(&LAPACKSupport::A, &mm, &nn, values, &mm,
	&wr[0], mu, &mm, mvt, &nn,
	&work[0], &lwork, &ipiv[0], &info);
  Assert (info==0, LAPACKSupport::ExcErrorCode("gesdd", info));
				   // Now resize work array and
  lwork = static_cast<int>(work[0] + .5);
#else
  int lwork = 3*std::min(mm,nn) +
	      std::max(std::max(mm,nn),4*std::min(mm,nn)*std::min(mm,nn)+4*std::min(mm,nn));
#endif
  work.resize(lwork);
				   // Do the actual SVD.
  gesdd(&LAPACKSupport::A, &mm, &nn, values, &mm,
	&wr[0], mu, &mm, mvt, &nn,
	&work[0], &lwork, &ipiv[0], &info);
  Assert (info==0, LAPACKSupport::ExcErrorCode("gesdd", info));

  work.resize(0);
  ipiv.resize(0);
  
  state = LAPACKSupport::svd;
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_inverse_svd(const double threshold)
{
  if (state == LAPACKSupport::matrix)
    compute_svd();
  
  Assert (state==LAPACKSupport::svd, ExcState(state));
  
  const double lim = wr[0]*threshold;
  for (unsigned int i=0;i<wr.size();++i)
    {
      if (wr[i] > lim)
	wr[i] = 1./wr[i];
      else
	wr[i] = 0.;
    }
  state = LAPACKSupport::inverse_svd;
}


template <typename number>
void
LAPACKFullMatrix<number>::invert()
{
  Assert(state == matrix || state == lu,
	 ExcState(state));
  const int mm = this->n_rows();
  const int nn = this->n_cols();
  Assert (nn == mm, ExcNotQuadratic());

  number* values = const_cast<number*> (this->data());
  ipiv.resize(mm);
  int info;

  if (state == matrix)
    {
      getrf(&mm, &nn, values, &mm, &ipiv[0], &info);

      Assert(info >= 0, ExcInternalError());
      Assert(info == 0, LACExceptions::ExcSingular());
    }

  inv_work.resize (mm);
  getri(&mm, values, &mm, &ipiv[0], &inv_work[0], &mm, &info);

  Assert(info >= 0, ExcInternalError());
  Assert(info == 0, LACExceptions::ExcSingular());

  state = inverse_matrix;
}


template <typename number>
void
LAPACKFullMatrix<number>::apply_lu_factorization(Vector<number>& v,
						 const bool transposed) const
{
  Assert(state == lu, ExcState(state));
  Assert(this->n_rows() == this->n_cols(),
	 LACExceptions::ExcNotQuadratic());

  const char* trans = transposed ? &T : &N;
  const int nn = this->n_cols();
  const number* values = this->data();
  int info;

  getrs(trans, &nn, &one, values, &nn, &ipiv[0],
	v.begin(), &nn, &info);

  Assert(info == 0, ExcInternalError());
}


template <typename number>
void
LAPACKFullMatrix<number>::compute_eigenvalues(
  const bool right,
  const bool left)
{
  Assert(state == matrix, ExcState(state));
  const int nn = this->n_cols();
  wr.resize(nn);
  wi.resize(nn);
  if (right) vr.resize(nn*nn);
  if (left)  vl.resize(nn*nn);

  number* values = const_cast<number*> (this->data());

  int info  = 0;
  int lwork = 1;
  const char * const jobvr = (right) ? (&V) : (&N);
  const char * const jobvl = (left)  ? (&V) : (&N);

				   // Optimal workspace query:

				   // The LAPACK routine DGEEV requires
				   // a sufficient large workspace variable,
				   // minimum requirement is
				   //    work.size>=4*nn.
				   // However, to improve performance, a
				   // somewhat larger workspace may be needed.

				   // SOME implementations of the LAPACK routine
				   // provide a workspace query call,
				   //   info:=0, lwork:=-1
				   // which returns an optimal value for the
				   // size of the workspace array
				   // (the PETSc 2.3.0 implementation does NOT
				   // provide this functionality!).

				   // define the DEAL_II_LIBLAPACK_NOQUERYMODE flag to
				   // disable the workspace query.
#ifndef DEAL_II_LIBLAPACK_NOQUERYMODE
  lwork = -1;
  work.resize(1);

  geev(jobvl, jobvr, &nn, values, &nn,
       &wr[0], &wi[0],
       &vl[0], &nn, &vr[0], &nn,
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
#else
  lwork = 4*nn;                    // no query mode
#endif
				   // resize workspace array
  work.resize((unsigned int) lwork);

				   // Finally compute the eigenvalues.
  geev(jobvl, jobvr, &nn, values, &nn,
       &wr[0], &wi[0],
       &vl[0], &nn, &vr[0], &nn,
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


template <typename number>
void
LAPACKFullMatrix<number>::compute_generalized_eigenvalues_symmetric (
  LAPACKFullMatrix<number> & B,
  std::vector<Vector<number> > & eigenvectors,
  const int itype)
{
  Assert(state == matrix, ExcState(state));
  const int nn = this->n_cols();
  Assert(static_cast<unsigned int>(nn) == this->n_rows(), ExcNotQuadratic());
  Assert(B.n_rows() == B.n_cols(), ExcNotQuadratic());
  Assert(static_cast<unsigned int>(nn) == B.n_cols(),
	 ExcDimensionMismatch (nn, B.n_cols()));
  Assert(eigenvectors.size() <= static_cast<unsigned int>(nn),
	 ExcMessage ("eigenvectors.size() > matrix.n_cols()"));

  wr.resize(nn);
  wi.resize(nn); //This is set purley for consistency reasons with the
                 //eigenvalues() function.

  number* values_A = const_cast<number*> (this->data());
  number* values_B = const_cast<number*> (B.data());

  int info  = 0;
  int lwork = 1;
  const char * const jobz = (eigenvectors.size() > 0) ? (&V) : (&N);
  const char * const uplo = (&U);

				   // Optimal workspace query:

				   // The LAPACK routine DSYGV requires
				   // a sufficient large workspace variable,
				   // minimum requirement is
				   //    work.size>=3*nn-1.
				   // However, to improve performance, a
				   // somewhat larger workspace may be needed.

				   // SOME implementations of the LAPACK routine
				   // provide a workspace query call,
				   //   info:=0, lwork:=-1
				   // which returns an optimal value for the
				   // size of the workspace array
				   // (the PETSc 2.3.0 implementation does NOT
				   // provide this functionality!).

				   // define the DEAL_II_LIBLAPACK_NOQUERYMODE flag to
				   // disable the workspace query.
#ifndef DEAL_II_LIBLAPACK_NOQUERYMODE
  lwork = -1;
  work.resize(1);

  sygv (&itype, jobz, uplo, &nn, values_A, &nn,
        values_B, &nn,
        &wr[0], &work[0], &lwork, &info);
				   // sygv returns info=0 on
				   // success. Since we only queried
				   // the optimal size for work,
				   // everything else would not be
				   // acceptable.
  Assert (info == 0, ExcInternalError());
				   // Allocate working array according
				   // to suggestion.
  lwork = (int) (work[0]+.1);
#else
  lwork = 3*nn-1 > 1 ? 3*nn-1 : 1; // no query mode
#endif
				   // resize workspace array
  work.resize((unsigned int) lwork);

				   // Finally compute the generalized
				   // eigenvalues.
  sygv (&itype, jobz, uplo, &nn, values_A, &nn,
        values_B, &nn,
        &wr[0], &work[0], &lwork, &info);
				   // Negative return value implies a
				   // wrong argument. This should be
				   // internal.

  Assert (info >=0, ExcInternalError());
  if (info != 0)
    std::cerr << "LAPACK error in sygv" << std::endl;

  for (unsigned int i=0; i < eigenvectors.size(); ++i)
  {
    unsigned int col_begin(i*nn);
    eigenvectors[i].reinit(nn, true);
    for (unsigned int j=0; j < static_cast<unsigned int>(nn); ++j)
    {
      eigenvectors[i](j) = values_A[col_begin+j];
    }
  }
  state = LAPACKSupport::State(eigenvalues | unusable);
}


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

template <typename number>
void
LAPACKFullMatrix<number>::print_formatted (
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
    } else {
      out.setf (std::ios::fixed, std::ios::floatfield);
      if (!width)
	width = precision+2;
    }

  for (unsigned int i=0; i<this->n_rows(); ++i)
    {
      for (unsigned int j=0; j<this->n_cols(); ++j)
	if (std::fabs(this->el(i,j)) > threshold)
	  out << std::setw(width)
	      << this->el(i,j) * denominator << ' ';
	else
	  out << std::setw(width) << zero_string << ' ';
      out << std::endl;
    };

  AssertThrow (out, ExcIO());
				   // reset output format
  out.flags (old_flags);
  out.precision(old_precision);
}


//----------------------------------------------------------------------//

template <typename number>
void
PreconditionLU<number>::initialize(const LAPACKFullMatrix<number>& M)
{
  matrix = &M;
  mem = 0;
}


template <typename number>
void
PreconditionLU<number>::initialize(const LAPACKFullMatrix<number>& M,
				    VectorMemory<Vector<number> >& V)
{
  matrix = &M;
  mem = &V;
}


template <typename number>
void
PreconditionLU<number>::vmult(Vector<number>& dst,
			       const Vector<number>& src) const
{
  dst = src;
  matrix->apply_lu_factorization(dst, false);
}


template <typename number>
void
PreconditionLU<number>::Tvmult(Vector<number>& dst,
				const Vector<number>& src) const
{
  dst = src;
  matrix->apply_lu_factorization(dst, true);
}


template <typename number>
void
PreconditionLU<number>::vmult(BlockVector<number>& dst,
			       const BlockVector<number>& src) const
{
  Assert(mem != 0, ExcNotInitialized());
  Vector<number>* aux = mem->alloc();
  *aux = src;
  matrix->apply_lu_factorization(*aux, false);
  dst = *aux;
}


template <typename number>
void
PreconditionLU<number>::Tvmult(BlockVector<number>& dst,
				const BlockVector<number>& src) const
{
  Assert(mem != 0, ExcNotInitialized());
  Vector<number>* aux = mem->alloc();
  *aux = src;
  matrix->apply_lu_factorization(*aux, true);
  dst = *aux;
}



template class LAPACKFullMatrix<double>;
template LAPACKFullMatrix<double> &
LAPACKFullMatrix<double>::operator = (const FullMatrix<double>& M);
template LAPACKFullMatrix<double> &
LAPACKFullMatrix<double>::operator = (const FullMatrix<float>& M);

template class LAPACKFullMatrix<float>;
template LAPACKFullMatrix<float> &
LAPACKFullMatrix<float>::operator = (const FullMatrix<double>& M);
template LAPACKFullMatrix<float> &
LAPACKFullMatrix<float>::operator = (const FullMatrix<float>& M);

template class LAPACKFullMatrix<long double>;
template LAPACKFullMatrix<long double> &
LAPACKFullMatrix<long double>::operator = (const FullMatrix<long double>& M);
template LAPACKFullMatrix<long double> &
LAPACKFullMatrix<long double>::operator = (const FullMatrix<double>& M);
template LAPACKFullMatrix<long double> &
LAPACKFullMatrix<long double>::operator = (const FullMatrix<float>& M);

template class PreconditionLU<double>;
template class PreconditionLU<float>;

DEAL_II_NAMESPACE_CLOSE
