// $Id$

// SparseMatrix template instantiation


/* Instantiation is controlled by preprocessor symbols:
 *
 * 1. TYPEMAT : numerical type used in the matrix
 * 2. TYPE2 : numerical type for function arguments
 */

#include <cmath>
#include <lac/sparse_matrix.templates.h>

#define TYPEMAT double

template class SparseMatrix<TYPEMAT>;

#define TYPE2 float

#include <lac/sparse_matrix.2.templates>

#undef TYPE2
#define TYPE2 double

#include <lac/sparse_matrix.2.templates>

#undef TYPE2
#define TYPE2 long double

#include <lac/sparse_matrix.2.templates>

#undef TYPE2
#undef TYPEMAT
