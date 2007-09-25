//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/filtered_matrix.templates.h>

DEAL_II_NAMESPACE_OPEN

template class FilteredMatrix<Vector<double> >;
template class FilteredMatrix<BlockVector<double> >;
template class FilteredMatrix<Vector<float> >;
template class FilteredMatrix<BlockVector<float> >;

DEAL_II_NAMESPACE_CLOSE
