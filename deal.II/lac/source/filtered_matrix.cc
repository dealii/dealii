//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/filtered_matrix.templates.h>

//TODO: Maybe split in several files?

#define TYPEMAT double
#define TYPEVEC double
#include "filtered_matrix.in.h"

template class FilteredMatrix<SparseMatrix<TYPEMAT>,Vector<TYPEVEC> >;
template class FilteredMatrix<BlockSparseMatrix<TYPEMAT>,BlockVector<TYPEVEC> >;

#undef TYPEVEC
#undef TYPEMAT



#define TYPEMAT float
#define TYPEVEC float
#include "filtered_matrix.in.h"

template class FilteredMatrix<SparseMatrix<TYPEMAT>,Vector<TYPEVEC> >;
template class FilteredMatrix<BlockSparseMatrix<TYPEMAT>,BlockVector<TYPEVEC> >;

#undef TYPEVEC
#undef TYPEMAT

