//---------------------------------------------------------------------------
//    $Id:  $
//    Version: $Name$
//
//    Copyright (C) 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// This file compiles the second half of the instantiations from fe_values.cc
// to get the memory consumption below 1.5gb with gcc.

#define FE_VALUES_INSTANTIATE_PART_TWO

#include "fe_values.cc"
