// $Id$

// SparseMatrix template instantiation


/* Instantiation is controlled by preprocessor symbols:
 *
 * 1. TYPEMAT : numerical type used in the matrix
 * 2. TYPE2 : numerical type for function arguments
 */

#include <cmath>
#include <lac/sparsematrix.templates.h>
#include <lac/sparse_vanka.templates.h>

#define TYPEMAT double

template class SparseMatrix<TYPEMAT>;
template class SparseVanka<TYPEMAT>;

#define TYPE2 float

#include <lac/sparsematrix.2.templates>

#undef TYPE2
#define TYPE2 double

#include <lac/sparsematrix.2.templates>

#undef TYPE2
#define TYPE2 long double

#include <lac/sparsematrix.2.templates>
