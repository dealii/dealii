//----------------------------------------------------------------
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
//----------------------------------------------------------------


// only compile this file if in 2d
#if deal_II_dimension == 2


#include <fe/fe_q.h>

// constraint matrices in 2d are now implemented by computing them on
// the fly for all polynomial degrees

template <>
const double * const 
FE_Q<2>::Matrices::constraint_matrices[] = {};


template <>
const unsigned int 
FE_Q<2>::Matrices::n_constraint_matrices = 0;



#else // #if deal_II_dimension
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // #if deal_II_dimension == 2
