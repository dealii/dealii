//----------------------------------------------------------------
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
//----------------------------------------------------------------



// only compile this file if in 1d
#if deal_II_dimension == 1


#include <fe/fe_q.h>

// Transfer matrices for finite elements

namespace FE_Q_1d
{
  static const double q1_into_q1_refined_0[] =
  {
	1., 0.,
	13.5/27., 13.5/27.,
  };

  static const double q1_into_q1_refined_1[] =
  {
	13.5/27., 13.5/27.,
	0., 1.,
  };

  static const double q2_into_q2_refined_0[] =
  {
	1., 0., 0.,
	0., 0., 1.,
	10.125/27., -3.375/27., 20.25/27.,
  };

  static const double q2_into_q2_refined_1[] =
  {
	0., 0., 1.,
	0., 1., 0.,
	-3.375/27., 10.125/27., 20.25/27.,
  };

  static const double q3_into_q3_refined_0[] =
  {
	1., 0., 0., 0.,
	-1.6875/27., -1.6875/27., 15.1875/27., 15.1875/27.,
	8.4375/27., 1.6875/27., 25.3125/27., -8.4375/27.,
	0., 0., 1., 0.
  };

  static const double q3_into_q3_refined_1[] =
  {
	-1.6875/27., -1.6875/27., 15.1875/27., 15.1875/27.,
	0., 1., 0., 0.,
	0., 0., 0., 1.,
	1.6875/27., 8.4375/27., -8.4375/27., 25.3125/27.,
  };

  static const double q4_into_q4_refined_0[] =
  {
	1., 0., 0., 0., 0.,
	0., 0., 0., 1., 0.,
	7.3828125/27., -1.0546875/27., 29.53125/27., -14.765625/27., 5.90625/27.,
	0., 0., 1., 0., 0.,
	-1.0546875/27., 0.6328125/27., 12.65625/27., 18.984375/27., -4.21875/27.,
  };

  static const double q4_into_q4_refined_1[] =
  {
	0., 0., 0., 1., 0.,
	0., 1., 0., 0., 0.,
	0.6328125/27., -1.0546875/27., -4.21875/27., 18.984375/27., 12.65625/27.,
	0., 0., 0., 0., 1.,
	-1.0546875/27., 7.3828125/27., 5.90625/27., -14.765625/27., 29.53125/27.,
  };
 
};   // namespace FE_Q_1d



// embedding matrices


const unsigned int  FE_Q<1>::Matrices::n_embedding_matrices = 4;

const double * const
FE_Q<1>::Matrices::embedding[][GeometryInfo<1>::children_per_cell] =
{
  {FE_Q_1d::q1_into_q1_refined_0, FE_Q_1d::q1_into_q1_refined_1},
  {FE_Q_1d::q2_into_q2_refined_0, FE_Q_1d::q2_into_q2_refined_1},
  {FE_Q_1d::q3_into_q3_refined_0, FE_Q_1d::q3_into_q3_refined_1},
  {FE_Q_1d::q4_into_q4_refined_0, FE_Q_1d::q4_into_q4_refined_1},
};


// No constraints in 1d
const unsigned int 
FE_Q<1>::Matrices::n_constraint_matrices = 0;

const double * const
FE_Q<1>::Matrices::constraint_matrices[] = { 0 };


#else // #if deal_II_dimension
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {}; };
#endif // #if deal_II_dimension == 1
