//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------



// only compile this file if in 1d. note that Raviart-Thomas elements
// do not make much sense in 1d anyway, so this file only contains
// dummy implementations to avoid linker errors due to missing symbols
#if deal_II_dimension == 1


#include <fe/fe_raviart_thomas.h>


template <>
const double * const
FE_RaviartThomas<1>::Matrices::embedding[][GeometryInfo<1>::children_per_cell] =
{};


template <>
const unsigned int
FE_RaviartThomas<1>::Matrices::n_embedding_matrices = 0;



// No constraints in 1d
template <>
const unsigned int 
FE_RaviartThomas<1>::Matrices::n_constraint_matrices = 0;


template <>
const double * const
FE_RaviartThomas<1>::Matrices::constraint_matrices[] = {};


#else // #if deal_II_dimension
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // #if deal_II_dimension == 1

