//------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//------------------------------------------------------------------------


#include <base/logstream.h>
#include <lac/sparse_matrix_ez.templates.h>
#include <iostream>

#define TYPEMAT double

template class SparseMatrixEZ<TYPEMAT>;
template void SparseMatrixEZ<TYPEMAT>::print_statistics(std::ostream&, bool);
template void SparseMatrixEZ<TYPEMAT>::print_statistics(LogStream&, bool);


#define TYPE2 float

#include <lac/sparse_matrix_ez.2.templates>

#undef TYPE2
#define TYPE2 double

#include <lac/sparse_matrix_ez.2.templates>

				 // a prerelease of gcc3.0 fails to
				 // compile this due to long double
//  #undef TYPE2
//  #define TYPE2 long double

//  #include <lac/sparse_matrix.2.templates>

#undef TYPE2
#undef TYPEMAT
