//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------



// only compile this file if in 1d. note that Nedelec elements do not
// make much sense in 1d anyway, so this file only contains dummy
// implementations to avoid linker errors due to missing symbols
#if deal_II_dimension == 1


#include <fe/fe_nedelec.h>

DEAL_II_NAMESPACE_OPEN

template <>
const double * const
FE_Nedelec<1>::Matrices::embedding[][GeometryInfo<1>::children_per_cell] =
{{0}};


template <>
const unsigned int
FE_Nedelec<1>::Matrices::n_embedding_matrices = 0;



// No constraints in 1d
template <>
const unsigned int 
FE_Nedelec<1>::Matrices::n_constraint_matrices = 0;


template <>
const double * const
FE_Nedelec<1>::Matrices::constraint_matrices[] = {0};

DEAL_II_NAMESPACE_CLOSE

#endif // #if deal_II_dimension == 1

