//----------------------------  sparse_matrix.float.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_matrix.float.cc  ---------------------------


#include <lac/sparse_matrix.templates.h>

#define TYPEMAT float

template class SparseMatrix<TYPEMAT>;

#define TYPEVEC float

#include "sparse_matrix.in.h"

#undef TYPEVEC
#define TYPEVEC double

#include "sparse_matrix.in.h"

				 // a prerelease of gcc3.0 fails to
				 // compile this due to long double
//  #undef TYPEVEC
//  #define TYPEVEC long double

//  #include <lac/sparse_matrix.2.templates>

#undef TYPEVEC
#undef TYPEMAT
