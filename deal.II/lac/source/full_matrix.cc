//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/full_matrix.templates.h>
#include <lac/lapack_templates.h>
#include <base/logstream.h>

DEAL_II_NAMESPACE_OPEN


                                   // Need to explicitly state the Lapack
                                   // inversion since it only works with
                                   // floats and doubles in case LAPACK was
                                   // detected by configure.
#if defined(HAVE_DGETRF_) && defined (HAVE_SGETRF_) && defined(HAVE_DGETRI_) && defined (HAVE_SGETRI_)

template <>
void
FullMatrix<float>::gauss_jordan ()
{
  Assert (!this->empty(), ExcEmptyMatrix());  
  Assert (this->n_cols() == this->n_rows(), ExcNotQuadratic());
  
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
  float* values = const_cast<float*> (this->data());
  ipiv.resize(nn);
  int info;

                                   // Use the LAPACK function getrf for 
                                   // calculating the LU factorization.
  getrf(&nn, &nn, values, &nn, &ipiv[0], &info);

  Assert(info >= 0, ExcInternalError());
  Assert(info == 0, LACExceptions::ExcSingular());

  inv_work.resize (nn);
                                   // Use the LAPACK function getri for
                                   // calculating the actual inverse using
                                   // the LU factorization.
  getri(&nn, values, &nn, &ipiv[0], &inv_work[0], &nn, &info);

  Assert(info >= 0, ExcInternalError());
  Assert(info == 0, LACExceptions::ExcSingular());
}

template <>
void
FullMatrix<double>::gauss_jordan ()
{
  Assert (!this->empty(), ExcEmptyMatrix());  
  Assert (this->n_cols() == this->n_rows(), ExcNotQuadratic());
  
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
  double* values = const_cast<double*> (this->data());
  ipiv.resize(nn);
  int info;

                                   // Use the LAPACK function getrf for 
                                   // calculating the LU factorization.
  getrf(&nn, &nn, values, &nn, &ipiv[0], &info);

  Assert(info >= 0, ExcInternalError());
  Assert(info == 0, LACExceptions::ExcSingular());

  inv_work.resize (nn);
                                   // Use the LAPACK function getri for
                                   // calculating the actual inverse using
                                   // the LU factorization.
  getri(&nn, values, &nn, &ipiv[0], &inv_work[0], &nn, &info);

  Assert(info >= 0, ExcInternalError());
  Assert(info == 0, LACExceptions::ExcSingular());
}

                                   // ... and now the usual instantiations
                                   // of gauss_jordan() and all the rest.
template void FullMatrix<long double>::gauss_jordan ();
template void FullMatrix<std::complex<float> >::gauss_jordan ();
template void FullMatrix<std::complex<double> >::gauss_jordan ();
template void FullMatrix<std::complex<long double> >::gauss_jordan ();

#else

template void FullMatrix<float>::gauss_jordan ();
template void FullMatrix<double>::gauss_jordan ();
template void FullMatrix<long double>::gauss_jordan ();
template void FullMatrix<std::complex<float> >::gauss_jordan ();
template void FullMatrix<std::complex<double> >::gauss_jordan ();
template void FullMatrix<std::complex<long double> >::gauss_jordan ();

#endif


#include "full_matrix.inst"


// do a few functions that currently don't fit the scheme because they have
// two template arguments that need to be different (the case of same
// arguments is covered by the default copy constructor and copy operator that
// is declared separately)

#define TEMPL_OP_EQ(S1,S2)			      \
  template FullMatrix<S1>& FullMatrix<S1>::operator = \
  (const FullMatrix<S2>&)

TEMPL_OP_EQ(double,float);
TEMPL_OP_EQ(float,double);

TEMPL_OP_EQ(long double,double);
TEMPL_OP_EQ(double,long double);

TEMPL_OP_EQ(long double,float);
TEMPL_OP_EQ(float,long double);


TEMPL_OP_EQ(std::complex<double>,std::complex<float>);
TEMPL_OP_EQ(std::complex<float>,std::complex<double>);

TEMPL_OP_EQ(std::complex<long double>,std::complex<double>);
TEMPL_OP_EQ(std::complex<double>,std::complex<long double>);

TEMPL_OP_EQ(std::complex<long double>,std::complex<float>);
TEMPL_OP_EQ(std::complex<float>,std::complex<long double>);

#undef TEMPL_OP_EQ


DEAL_II_NAMESPACE_CLOSE
