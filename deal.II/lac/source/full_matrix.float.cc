// $Id$

// Driver for FullMatrix template instantiation.

/* Instantiation is controlled by preprocessor symbols:
 *
 * 1. TYPEMAT : numerical type used in the matrix
 * 2. TYPEVEC : numerical type for vector entries
 * 3. TYPERES : numerical type for entries in the right hand side vector
 * 4. TYPEMAT2: numerical type for the second matrix
 */

#include <lac/fullmatrix.templates.h>

#define TYPEMAT float

template class FullMatrix<TYPEMAT>;

#define TYPEMAT2 float

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

template void FullMatrix<TYPEMAT>::vmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template void FullMatrix<TYPEMAT>::Tvmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template double FullMatrix<TYPEMAT>::residual(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const Vector<TYPERES>&) const;
template double FullMatrix<TYPEMAT>::matrix_norm (const Vector<TYPEVEC> &) const;
template double FullMatrix<TYPEMAT>::matrix_scalar_product(const Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::forward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::backward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::householder(Vector<TYPEVEC>&);
template double FullMatrix<TYPEMAT>::least_squares(Vector<TYPEVEC>&, Vector<TYPEVEC>&);
template void FullMatrix<TYPEMAT>::gsmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const iVector&) const;

#undef TYPEVEC
#define TYPEVEC float

template void FullMatrix<TYPEMAT>::vmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template void FullMatrix<TYPEMAT>::Tvmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const bool) const;
template double FullMatrix<TYPEMAT>::residual(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const Vector<TYPERES>&) const;
template double FullMatrix<TYPEMAT>::matrix_norm (const Vector<TYPEVEC> &) const;
template double FullMatrix<TYPEMAT>::matrix_scalar_product(const Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::forward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::backward(Vector<TYPEVEC>&, const Vector<TYPEVEC>&) const;
template void FullMatrix<TYPEMAT>::householder(Vector<TYPEVEC>&);
template double FullMatrix<TYPEMAT>::least_squares(Vector<TYPEVEC>&, Vector<TYPEVEC>&);
template void FullMatrix<TYPEMAT>::gsmult(Vector<TYPEVEC>&, const Vector<TYPEVEC>&, const iVector&) const;

