//----------------------------  sparse_matrix.float.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_matrix.float.cc  ---------------------------


#include <cmath>
#include <lac/sparse_matrix.templates.h>


#define TYPEMAT float

template class SparseMatrix<TYPEMAT>;


#define TYPE2 float

#include <lac/sparse_matrix.2.templates>

#undef TYPE2
#define TYPE2 double

#include <lac/sparse_matrix.2.templates>

#undef TYPE2
#define TYPE2 long double

#include <lac/sparse_matrix.2.templates>
