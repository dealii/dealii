// $Id$

// SparseMatrix template instantiation


/* Instantiation is controlled by preprocessor symbols:
 *
 * 1. TYPEMAT : numerical type used in the matrix
 * 2. TYPE2 : numerical type for function arguments
 */

#include <cmath>
#include <lac/sparsematrix.templates.h>

#define TYPEMAT float

template class SparseMatrix<TYPEMAT>;

#define TYPE2 float

template SparseMatrix<TYPEMAT> &
SparseMatrix<TYPEMAT>::copy_from (const SparseMatrix<TYPE2> &);

template void SparseMatrix<TYPEMAT>::add_scaled (const TYPEMAT,
						 const SparseMatrix<TYPE2> &);

template void SparseMatrix<TYPEMAT>::vmult (Vector<TYPE2> &,
					    const Vector<TYPE2> &) const;
template void SparseMatrix<TYPEMAT>::Tvmult (Vector<TYPE2> &,
					     const Vector<TYPE2> &) const;

template double
SparseMatrix<TYPEMAT>::matrix_norm (const Vector<TYPE2> &) const;

template double SparseMatrix<TYPEMAT>::residual (Vector<TYPE2> &,
					       const Vector<TYPE2> &,
					       const Vector<TYPE2> &) const;

template void SparseMatrix<TYPEMAT>::precondition_SSOR (Vector<TYPE2> &,
						      const Vector<TYPE2> &,
						      TYPEMAT) const;

template void SparseMatrix<TYPEMAT>::precondition_SOR (Vector<TYPE2> &,
						     const Vector<TYPE2> &,
						     TYPEMAT) const;

template void SparseMatrix<TYPEMAT>::precondition_Jacobi (Vector<TYPE2> &,
							const Vector<TYPE2> &,
							TYPEMAT) const;

template void SparseMatrix<TYPEMAT>::SOR (Vector<TYPE2> &, TYPEMAT) const;
template void SparseMatrix<TYPEMAT>::SSOR (Vector<TYPE2> &, TYPEMAT) const;

#undef TYPE2
#define TYPE2 double

template SparseMatrix<TYPEMAT> &
SparseMatrix<TYPEMAT>::copy_from (const SparseMatrix<TYPE2> &);

template void SparseMatrix<TYPEMAT>::add_scaled (const TYPEMAT,
						 const SparseMatrix<TYPE2> &);

template void SparseMatrix<TYPEMAT>::vmult (Vector<TYPE2> &,
					    const Vector<TYPE2> &) const;
template void SparseMatrix<TYPEMAT>::Tvmult (Vector<TYPE2> &,
					     const Vector<TYPE2> &) const;

template double
SparseMatrix<TYPEMAT>::matrix_norm (const Vector<TYPE2> &) const;

template double SparseMatrix<TYPEMAT>::residual (Vector<TYPE2> &,
					       const Vector<TYPE2> &,
					       const Vector<TYPE2> &) const;

template void SparseMatrix<TYPEMAT>::precondition_SSOR (Vector<TYPE2> &,
						      const Vector<TYPE2> &,
						      TYPEMAT) const;

template void SparseMatrix<TYPEMAT>::precondition_SOR (Vector<TYPE2> &,
						     const Vector<TYPE2> &,
						     TYPEMAT) const;

template void SparseMatrix<TYPEMAT>::precondition_Jacobi (Vector<TYPE2> &,
							const Vector<TYPE2> &,
							TYPEMAT) const;

template void SparseMatrix<TYPEMAT>::SOR (Vector<TYPE2> &, TYPEMAT) const;
template void SparseMatrix<TYPEMAT>::SSOR (Vector<TYPE2> &, TYPEMAT) const;
