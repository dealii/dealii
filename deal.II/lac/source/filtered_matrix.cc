//----------------------------  filtered_matrix.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  filtered_matrix.cc  ---------------------------


#include <lac/filtered_matrix.templates.h>


template class FilteredMatrix<SparseMatrix<double>,Vector<double> >;
template class FilteredMatrix<BlockSparseMatrix<double>,BlockVector<double> >;
