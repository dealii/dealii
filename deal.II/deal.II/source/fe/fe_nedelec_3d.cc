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

// Transfer matrices for finite elements


// only compile this file if in 3d
#if deal_II_dimension == 3

#include <fe/fe_nedelec.h>

// Transfer matrices for finite elements: have one matrix for each of
// the four child cells which tells us how the degrees of freedom on
// the child cell are obtained from the degrees of freedom on the
// mother cell
//
// TODO: [Anna] check whether the following paragraph is correct. if so, then please multiply the values in the eight following matrices by two

// note the following: since the shape functions themselves and not
// only the gradients are transformed using the mapping object from
// the unit cell to the real cell, the actual values of the function
// on the real cell is degree of freedom times value of the shape
// function on the unit cell times Jacobian. Thus, what has the DoF
// value 1 on the mother cell must have the DoF value 2 on the child
// cell since the latter is smaller by a (linear scaling) factor of
// two.
namespace FE_Nedelec_3d
{
  static const double q1_into_q1_refined_0[] =
  {
	1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	0., 0.5, 0., 0.5, 0., 0., 0., 0.,0.,0.,0.,0.,
	0.5, 0., 0.5, 0., 0., 0., 0., 0.,0.,0.,0.,0.,
	0., 0., 0., 1., 0., 0., 0., 0.,0.,0.,0.,0.,
	0.5, 0., 0., 0., 0.5, 0., 0., 0.,0.,0.,0.,0.,
	0., 0.25, 0., 0.25, 0., 0.25, 0., 0.25,0.,0.,0.,0.,
	0.25, 0., 0.25, 0., 0.25, 0., 0.25, 0.,0.,0.,0.,0.,
	0., 0., 0., 0.5, 0., 0., 0., 0.5,0.,0.,0.,0.,
	0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5, 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0.25, 0.25, 0.25,
	0., 0., 0., 0., 0., 0., 0., 0., 0.5, .0, 0., 0.5,
  };

  static const double q1_into_q1_refined_1[] =
  {
	1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	0.5, 0., 0.5, 0., 0., 0., 0., 0.,0.,0.,0.,0.,
	0., 0.5, 0., 0.5, 0., 0., 0., 0.,0.,0.,0.,0.,
	0.5, 0., 0., 0., 0.5, 0., 0., 0.,0.,0.,0.,0.,
	0., 0.5, 0., 0., 0., 0.5, 0., 0.,0.,0.,0.,0.,
	0.25, 0., 0.25, 0., 0.25, 0., 0.25, 0.,0.,0.,0.,0.,
	0., 0.25, 0., 0.25, 0., 0.25, 0., 0.25,0.,0.,0.,0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5, 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5, 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0.25, 0.25, 0.25,
  };

  static const double q1_into_q1_refined_2[] =
  {
	0.5, 0., 0.5, 0., 0., 0., 0., 0.,0.,0.,0.,0.,
	0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	0., 0.5, 0., 0.5, 0., 0., 0., 0.,0.,0.,0.,0.,
	0.25, 0., 0.25, 0., 0.25, 0., 0.25, 0.,0.,0.,0.,0.,
	0., 0.5, 0., 0., 0., 0.5, 0., 0.,0.,0.,0.,0.,
	0., 0., 0.5, 0., 0., 0., 0.5, 0.,0.,0.,0.,0.,
	0., 0.25, 0., 0.25, 0., 0.25, 0., 0.25,0.,0.,0.,0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0.25, 0.25, 0.25,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5, 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5,
  };

  static const double q1_into_q1_refined_3[] =
  {
	0.5, 0., 0.5, 0., 0., 0., 0., 0.,0.,0.,0.,0.,
	0., 0.5, 0., 0.5, 0., 0., 0., 0.,0.,0.,0.,0.,
	0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
	0.25, 0., 0.25, 0., 0.25, 0., 0.25, 0.,0.,0.,0.,0.,
	0., 0.25, 0., 0.25, 0., 0.25, 0., 0.25,0.,0.,0.,0.,
	0., 0., 0.5, 0., 0., 0., 0.5, 0.,0.,0.,0.,0.,
	0., 0., 0., 0.5, 0., 0., 0., 0.5,0.,0.,0.,0.,
	0., 0., 0., 0., 0., 0., 0., 0.,0.5,0.,0.,0.5,
	0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0.25, 0.25, 0.25,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
  };

  static const double q1_into_q1_refined_4[] =
  {
	0.5, 0., 0., 0., 0.5, 0., 0., 0.,0.,0.,0.,0.,
	0., 0.25, 0., 0.25, 0., 0.25, 0., 0.25,0.,0.,0.,0.,
	0.25, 0., 0.25, 0., 0.25, 0., 0.25, 0.,0.,0.,0.,0.,
	0., 0., 0., 0.5, 0., 0., 0., 0.5, 0., 0., 0., 0.,
	0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 0.5, 0., 0.5,0.,0.,0.,0.,
	0., 0., 0., 0., 0.5, 0., 0.5, 0.,0.,0.,0.,0.,
	0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5, 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0.25, 0.25, 0.25,
	0., 0., 0., 0., 0., 0., 0., 0.,0.5,0.,0.,0.5,
  };

  static const double q1_into_q1_refined_5[] =
  { 
	0.5, 0., 0., 0., 0.5, 0., 0., 0.,0.,0.,0.,0.,
	0., 0.5, 0., 0., 0., 0.5, 0., 0.,0.,0.,0.,0.,
	0.25, 0., 0.25, 0., 0.25, 0., 0.25, 0.,0.,0.,0.,0.,
	0., 0.25, 0., 0.25, 0., 0.25, 0., 0.25,0.,0.,0.,0.,
	0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0.5, 0., 0.5, 0.,0.,0.,0.,0.,
	0., 0., 0., 0., 0., 0.5, 0., 0.5,0.,0.,0.,0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5, 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5, 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0.25, 0.25, 0.25,

  };

  static const double q1_into_q1_refined_6[] =
  {
	0.25, 0., 0.25, 0., 0.25, 0., 0.25, 0.,0.,0.,0.,0.,
	0., 0.5, 0., 0., 0., 0.5, 0., 0.,0.,0.,0.,0.,
	0., 0., 0.5, 0., 0., 0., 0.5, 0.,0.,0.,0.,0.,
	0., 0.25, 0., 0.25, 0., 0.25, 0., 0.25,0.,0.,0.,0.,
	0., 0., 0., 0., 0.5, 0., 0.5, 0.,0.,0.,0.,0.,
	0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 0.5, 0., 0.5, 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0.25, 0.25, 0.25,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5, 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5,
	

  };

  static const double q1_into_q1_refined_7[] =
  {
	0.25, 0., 0.25, 0., 0.25, 0., 0.25, 0.,0.,0.,0.,0.,
	0., 0.25, 0., 0.25, 0., 0.25, 0., 0.25,0.,0.,0.,0.,
	0., 0., 0.5, 0., 0., 0., 0.5, 0.,0.,0.,0.,0.,
	0., 0., 0., 0.5, 0., 0., 0., 0.5,0.,0.,0.,0.,
	0., 0., 0., 0., 0.5, 0., 0.5, 0.,0.,0.,0.,0.,
	0., 0., 0., 0., 0., 0.5, 0., 0.5,0.,0.,0.,0.,
	0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0., 0., 0.5,
	0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0.25, 0.25, 0.25,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5,
	0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.1,
  };

};  // namespace FE_Nedelec_3d


// embedding matrices

template <>
const double * const 
FE_Nedelec<3>::Matrices::embedding[][GeometryInfo<3>::children_per_cell] =
{
      { FE_Nedelec_3d::q1_into_q1_refined_0, FE_Nedelec_3d::q1_into_q1_refined_1,
	FE_Nedelec_3d::q1_into_q1_refined_2, FE_Nedelec_3d::q1_into_q1_refined_3,
	FE_Nedelec_3d::q1_into_q1_refined_4, FE_Nedelec_3d::q1_into_q1_refined_5,
	FE_Nedelec_3d::q1_into_q1_refined_6, FE_Nedelec_3d::q1_into_q1_refined_7 }
};


template <>
const unsigned int
FE_Nedelec<3>::Matrices::n_embedding_matrices
= sizeof(FE_Nedelec<3>::Matrices::embedding) /
sizeof(FE_Nedelec<3>::Matrices::embedding[0]);



// Constraint matrices: how do the new value on child faces depend on
// the values on the mother face if that face has a hanging node
//
// Here, the same applies as for the embedding matrices: since the DoF
// values are not only multiplied by the values of the shape function
// on the unit cell, but also by the transformation, we have to
// multiply the value on the large face by 1/2 to get the same value
// back on the small face
namespace FE_Nedelec_3d 
{
  static const double constraint_q1[] =
  {
	0, .25, 0, .25,  // first the four interior lines
	.25, 0, .25, 0,
	0, .25, 0, .25,
	.25, 0, .25, 0,
	.5, 0, 0, 0,  // then the two child lines of each of the four outer
	.5, 0, 0, 0,  // ones. since the shape functions are constant on each
	0, .5, 0, 0,  // line, the two child lines get the same weights, modulo
	0, .5, 0, 0,  // the issue with the division by length scaling
	0, 0, .5, 0,
	0, 0, .5, 0,
	0, 0, 0, .5,
	0, 0, 0, .5
  };
};



template <>
const double * const 
FE_Nedelec<3>::Matrices::constraint_matrices[] =
{
      FE_Nedelec_3d::constraint_q1
};



template <>
const unsigned int 
FE_Nedelec<3>::Matrices::n_constraint_matrices
= sizeof(FE_Nedelec<3>::Matrices::constraint_matrices) /
sizeof(FE_Nedelec<3>::Matrices::constraint_matrices[0]);



#else // #if deal_II_dimension
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {}; };
#endif // #if deal_II_dimension == 3
