//----------------------------  sparse_matrix.double.cc  ---------------------------
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
//----------------------------  sparse_matrix.double.cc  ---------------------------


#include <cmath>
#include <lac/sparse_matrix.templates.h>

#define TYPEMAT double

template class SparseMatrix<TYPEMAT>;

#define TYPE2 float

#include <lac/sparse_matrix.2.templates>

#undef TYPE2
#define TYPE2 double

#include <lac/sparse_matrix.2.templates>

				 // a prerelease of gcc3.0 fails to
				 // compile this due to long double
//  #undef TYPE2
//  #define TYPE2 long double

//  #include <lac/sparse_matrix.2.templates>

#undef TYPE2
#undef TYPEMAT
