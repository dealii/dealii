//----------------------------  full_matrix.double.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix.double.cc  ---------------------------


#include <lac/full_matrix.templates.h>
#include <base/logstream.h>

#define TYPEMAT double

template class FullMatrix<TYPEMAT>;

#define TYPEMAT2 double

//template FullMatrix<TYPEMAT>& FullMatrix<TYPEMAT>::operator =(const FullMatrix<TYPEMAT2>&);
template void FullMatrix<TYPEMAT>::fill (const FullMatrix<TYPEMAT2>&, const unsigned, const unsigned);
template void FullMatrix<TYPEMAT>::reinit (const FullMatrix<TYPEMAT2>&);
template void FullMatrix<TYPEMAT>::add (const TYPEMAT, const FullMatrix<TYPEMAT2>&);
template void FullMatrix<TYPEMAT>::Tadd (const TYPEMAT, const FullMatrix<TYPEMAT2>&);
template void FullMatrix<TYPEMAT>::mmult (FullMatrix<TYPEMAT2>&, const FullMatrix<TYPEMAT2>&) const;
template void FullMatrix<TYPEMAT>::Tmmult (FullMatrix<TYPEMAT2>&, const FullMatrix<TYPEMAT2>&) const;
template void FullMatrix<TYPEMAT>::add_diag (const TYPEMAT, const FullMatrix<TYPEMAT2>&);


#define TYPEVEC double
#define TYPERES double

template void FullMatrix<TYPEMAT>::fill (const TYPEVEC*);
template void FullMatrix<TYPEMAT>::vmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template void FullMatrix<TYPEMAT>::Tvmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template double FullMatrix<TYPEMAT>::residual(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const Vector<TYPERES>&) const;
template TYPEVEC FullMatrix<TYPEMAT>::matrix_norm_square (const Vector<TYPEVEC> &) const;
template TYPEVEC FullMatrix<TYPEMAT>::matrix_scalar_product(const Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::forward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::backward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::householder(Vector<TYPEVEC>&);
template double FullMatrix<TYPEMAT>::least_squares(Vector<TYPEVEC>&, Vector<TYPEVEC>&);

#undef TYPEVEC
#define TYPEVEC float

template void FullMatrix<TYPEMAT>::fill (const TYPEVEC*);
template void FullMatrix<TYPEMAT>::vmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template void FullMatrix<TYPEMAT>::Tvmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template double FullMatrix<TYPEMAT>::residual(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const Vector<TYPERES>&) const;
template TYPEVEC FullMatrix<TYPEMAT>::matrix_norm_square (const Vector<TYPEVEC> &) const;
template TYPEVEC FullMatrix<TYPEMAT>::matrix_scalar_product(const Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::forward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::backward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::householder(Vector<TYPEVEC>&);
template double FullMatrix<TYPEMAT>::least_squares(Vector<TYPEVEC>&, Vector<TYPEVEC>&);

#undef TYPERES
#define TYPERES float

template double FullMatrix<TYPEMAT>::residual(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const Vector<TYPERES>&) const;

// Experimental code

#ifdef HAVE_LIBLAPACK
extern "C" int dgels_ (const char* trans,
		       const unsigned int* M, const unsigned int* N,
		       const unsigned int* NRHS,
		       double* A, const unsigned int* LDA,
		       double* B, const unsigned int* LDB,
		       double* WORK, const unsigned int* LWORK,
		       int* INFO);
extern "C" int dgelss_ (const unsigned int* M, const unsigned int* N,
			const unsigned int* NRHS,
			double* A, const unsigned int* LDA,
			double* B, const unsigned int* LDB,
			double* S, const double* RCOND,
			int* RANK,
			double* WORK, const unsigned int* LWORK,
			int* INFO);


template<>
void
FullMatrix<double>::invert (const FullMatrix<double> &M)
{
  Assert (val != 0, ExcEmptyMatrix());
  
  Assert (dim_range == dim_image, ExcNotQuadratic());
  Assert (dim_range == M.dim_range,
          ExcDimensionMismatch(dim_range,M.dim_range));
  Assert (dim_image == M.dim_image,
	  ExcDimensionMismatch(dim_image,M.dim_image));

  clear();
  diagadd(1.);

  const unsigned int lwork = 10*dim_range*dim_range;
  double* work = new double[lwork];
  int info;
  
//  const char* trans = "N";
  int rank;
  const double rcond=-1.;
  double* s = new double[dim_range];
  
  double* matrix = new double[dim_range*dim_range];
  std::copy (&M.val[0], &M.val[dim_image*dim_range], matrix);
  

  int erg = dgelss_ (&dim_range, &dim_range, &dim_range,
		     matrix, &dim_range,
		     val, &dim_range,
		     s, &rcond, &rank,
		     work, &lwork,
		     &info);
//    int erg = dgels_ (trans, &dim_range, &dim_range, &dim_range,
//  		    M.val, &dim_range,
//  		    val, &dim_range,
//  		    work, &lwork,
//  		    &info);

//  double condition = s[0]/s[dim_range-1];
  
  if (info!=0)
    deallog << "inverting error " << info << ' ' << erg << endl;
  if (rank<(int)dim_range)
    deallog << "rank deficiency " << rank << endl;
  delete[] work;
  delete[] s;
  delete[] matrix;
}

#endif
