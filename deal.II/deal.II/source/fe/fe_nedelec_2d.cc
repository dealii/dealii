//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------


// only compile this file if in 2d
#if deal_II_dimension == 2


#include <fe/fe_nedelec.h>

// Transfer matrices for finite elements: have one matrix for each of
// the four child cells which tells us how the degrees of freedom on
// the child cell are obtained from the degrees of freedom on the
// mother cell
namespace FE_Nedelec_2d
{
  static const double q1_into_q1_refined_0[] =
  {
	1.,  0,  0,  0,
	0,   0.5,0,  0.5,
	0.5, 0,  0.5,0,
	0,   0,  0,  1 
  };

  static const double q1_into_q1_refined_1[] =
  {
//[Anna]:done	
  	1., 0., 0., 0.,
  	0., 1., 0., 0.,
  	0.5, 0., 0.5, 0.,
  	 0., 0.5, 0., 0.5,
  };

  static const double q1_into_q1_refined_2[] =
  {
  	0.5, 0., 0.5, 0.,
 	0., 1., 0., 0.,
	0., 0., 1., 0.,
  	0., 0.5, 0., 0.5,
  };

  static const double q1_into_q1_refined_3[] =
  {
  	0.5, 0., 0.5, 0.,
  	0., 0.5, 0., 0.5,
	0., 0., 1., 0.,
  	0., 0., 0., 1.,
  };
};  // namespace FE_Nedelec_2d


// embedding matrices

template <>
const double * const 
FE_Nedelec<2>::Matrices::embedding[][GeometryInfo<2>::children_per_cell] =
{
  { FE_Nedelec_2d::q1_into_q1_refined_0, FE_Nedelec_2d::q1_into_q1_refined_1,
      FE_Nedelec_2d::q1_into_q1_refined_2, FE_Nedelec_2d::q1_into_q1_refined_3 }
};


template <>
const unsigned int
FE_Nedelec<2>::Matrices::n_embedding_matrices
  = sizeof(FE_Nedelec<2>::Matrices::embedding) /
    sizeof(FE_Nedelec<2>::Matrices::embedding[0]);


// Constraint matrices: how do the new value on child faces depend on
// the values on the mother face if that face has a hanging node
namespace FE_Nedelec_2d 
{
  static const double constraint_q1[] =
  {
					 // the function is constant
					 // along each edge, so each
					 // degree of freedom on the
					 // refined edge has the same
					 // value as that on the
					 // coarse edge
  	1., 1.     
  };

};


template <>
const double * const 
FE_Nedelec<2>::Matrices::constraint_matrices[] =
{
  FE_Nedelec_2d::constraint_q1
};


template <>
const unsigned int 
FE_Nedelec<2>::Matrices::n_constraint_matrices
  = sizeof(FE_Nedelec<2>::Matrices::constraint_matrices) /
    sizeof(FE_Nedelec<2>::Matrices::constraint_matrices[0]);



#else // #if deal_II_dimension
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {}; };
#endif // #if deal_II_dimension == 2
