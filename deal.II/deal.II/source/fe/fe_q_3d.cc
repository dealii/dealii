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

// Transfer matrices for finite elements


// only compile this file if in 3d
#if deal_II_dimension == 3

#include <fe/fe_q.h>

// Constraint matrices taken from Wolfgangs old version

namespace FE_Q_3d 
{
  static const double constraint_q1[] =
  {
	.25,.25,.25,.25,
	.5,.5,0.,0.,
	0.,.5,.5,0.,
	0.,0.,.5,.5,
	.5,0.,0.,.5
  };

  static const double constraint_q2[] =
  {
	0,0,0,0,0,0,0,0,1,
	0,0,0,0,1,0,0,0,0,
	0,0,0,0,0,1,0,0,0,
	0,0,0,0,0,0,1,0,0,
	0,0,0,0,0,0,0,1,0,
	0,0,0,0,.375,0,-.125,0,.75,
	0,0,0,0,0,.375,0,-.125,.75,
	0,0,0,0,-.125,0,.375,0,.75,
	0,0,0,0,0,-.125,0,.375,.75,
	.375,-.125,0,0,.75,0,0,0,0,
	-.125,.375,0,0,.75,0,0,0,0,
	0,.375,-.125,0,0,.75,0,0,0,
	0,-.125,.375,0,0,.75,0,0,0,
	0,0,-.125,.375,0,0,.75,0,0,
	0,0,.375,-.125,0,0,.75,0,0,
	.375,0,0,-.125,0,0,0,.75,0,
	-.125,0,0,.375,0,0,0,.75,0,
	.140625,-.046875,.015625,-.046875,.28125,-.09375,-.09375,.28125,.5625,
	-.046875,.140625,-.046875,.015625,.28125,.28125,-.09375,-.09375,.5625,
	.015625,-.046875,.140625,-.046875,-.09375,.28125,.28125,-.09375,.5625,
	-.046875,.015625,-.046875,.140625,-.09375,-.09375,.28125,.28125,.5625
  };
}



template <>
const double * const 
FE_Q<3>::Matrices::constraint_matrices[] =
{
  FE_Q_3d::constraint_q1,
  FE_Q_3d::constraint_q2
};



template <>
const unsigned int 
FE_Q<3>::Matrices::n_constraint_matrices
  = sizeof(FE_Q<3>::Matrices::constraint_matrices) /
    sizeof(FE_Q<3>::Matrices::constraint_matrices[0]);



#else // #if deal_II_dimension
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // #if deal_II_dimension == 3
