//----------------------------  full_matrix.float.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix.float.cc  ---------------------------


#include <lac/full_matrix.templates.h>

#define TYPEMAT float

template class FullMatrix<TYPEMAT>;

#define TYPEMAT2 float

//template FullMatrix<TYPEMAT>& FullMatrix<TYPEMAT>::operator =(const FullMatrix<TYPEMAT2>&);
template void FullMatrix<TYPEMAT>::fill (const FullMatrix<TYPEMAT2>&, const unsigned, const unsigned);
template void FullMatrix<TYPEMAT>::add (const TYPEMAT, const FullMatrix<TYPEMAT2>&);
template void FullMatrix<TYPEMAT>::Tadd (const TYPEMAT, const FullMatrix<TYPEMAT2>&);
template void FullMatrix<TYPEMAT>::mmult (FullMatrix<TYPEMAT2>&, const FullMatrix<TYPEMAT2>&, const bool) const;
template void FullMatrix<TYPEMAT>::Tmmult (FullMatrix<TYPEMAT2>&, const FullMatrix<TYPEMAT2>&, const bool) const;
template void FullMatrix<TYPEMAT>::add_diag (const TYPEMAT, const FullMatrix<TYPEMAT2>&);


#define TYPEVEC double
#define TYPERES double

template void FullMatrix<TYPEMAT>::fill_permutation (const FullMatrix<TYPEVEC>&,
						     const std::vector<unsigned int>&,
						     const std::vector<unsigned int>&);
template void FullMatrix<TYPEMAT>::vmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template void FullMatrix<TYPEMAT>::Tvmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template double FullMatrix<TYPEMAT>::residual(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const Vector<TYPERES>&) const;
template TYPEVEC FullMatrix<TYPEMAT>::matrix_norm_square (const Vector<TYPEVEC> &) const;
template TYPEVEC FullMatrix<TYPEMAT>::matrix_scalar_product(const Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::forward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::backward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::householder(Vector<TYPEVEC>&);
template double FullMatrix<TYPEMAT>::least_squares(Vector<TYPEVEC>&, Vector<TYPEVEC>&);
template
void FullMatrix<TYPEMAT>::precondition_Jacobi (Vector<TYPEVEC> &,
					       const Vector<TYPEVEC> &,
					       const TYPEMAT) const;

#undef TYPEVEC
#define TYPEVEC float

template void FullMatrix<TYPEMAT>::fill_permutation (const FullMatrix<TYPEVEC>&,
						     const std::vector<unsigned int>&,
						     const std::vector<unsigned int>&);
template void FullMatrix<TYPEMAT>::vmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template void FullMatrix<TYPEMAT>::Tvmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template double FullMatrix<TYPEMAT>::residual(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const Vector<TYPERES>&) const;
template TYPEVEC FullMatrix<TYPEMAT>::matrix_norm_square (const Vector<TYPEVEC> &) const;
template TYPEVEC FullMatrix<TYPEMAT>::matrix_scalar_product(const Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::forward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::backward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::householder(Vector<TYPEVEC>&);
template double FullMatrix<TYPEMAT>::least_squares(Vector<TYPEVEC>&, Vector<TYPEVEC>&);
template
void FullMatrix<TYPEMAT>::precondition_Jacobi (Vector<TYPEVEC> &,
					       const Vector<TYPEVEC> &,
					       const TYPEMAT) const;


#undef TYPERES
#define TYPERES float

template double FullMatrix<TYPEMAT>::residual(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const Vector<TYPERES>&) const;
