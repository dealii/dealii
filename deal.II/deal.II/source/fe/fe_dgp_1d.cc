//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------


// only compile in 1d
#if deal_II_dimension == 1

#include <fe/fe_dgp.h>

// Transfer matrices for finite elements

#define SQRT3 1.732050807569

namespace FE_DGP_1d
{
  static const double dgp0_into_dgp0_refined_0[] =
  {
	1.,
  };

  static const double dgp0_into_dgp0_refined_1[] =
  {
	1.,
  };

  static const double dgp1_into_dgp1_refined_0[] =
  {
	1., -SQRT3/2.,
	0, .5,
  };

  static const double dgp1_into_dgp1_refined_1[] =
  {
	1., SQRT3/2.,
	0, .5,
  };

  static const double dgp2_into_dgp2_refined_0[] =
  {
	1., -SQRT3/2., 0,
	0, .5, -26.14263759/27.,
	0, 0, .25,
  };

  static const double dgp2_into_dgp2_refined_1[] =
  {
	1., SQRT3/2., 0,
	0, .5, 26.14263759/27.,
	0, 0, .25,
  };

  static const double dgp3_into_dgp3_refined_0[] =
  {
	1., -SQRT3/2., 0, 8.929410675/27.,
	0, .5, -26.14263759/27., 15.46619297/27.,
	0, 0, .25, -19.96676927/27.,
	0, 0, 0, .125,
  };

  static const double dgp3_into_dgp3_refined_1[] =
  {
	1., SQRT3/2., 0, -8.929410675/27.,
	0, .5, 26.14263759/27., 15.46619297/27.,
	0, 0, .25, 19.96676927/27.,
	0, 0, 0, .125,
  };

  static const double dgp4_into_dgp4_refined_0[] =
  {
	1., -SQRT3/2., 0, 8.929410675/27., 0,
	0, .5, -26.14263759/27., 15.46619297/27., SQRT3/8.,
	0, 0, .25, -19.96676927/27., 22.64018827/27.,
	0, 0, 0, .125, -13.39411601/27.,
	0, 0, 0, 0, 1./16.,
  };

  static const double dgp4_into_dgp4_refined_1[] =
  {
	1., SQRT3/2., 0, -8.929410675/27., 0.,
	0, .5, 26.14263759/27., 15.46619297/27., -SQRT3/8.,
	0, 0, .25, 19.96676927/27., 22.64018827/27.,
	0, 0, 0, .125, 13.39411601/27.,
	0, 0, 0, 0, 1./16.,
  };

  static const double dgp5_into_dgp5_refined_0[] =
  {
	1., -SQRT3/2., 0, 8.929410675/27., 0, -5.596804334/27.,
	0, .5, -26.14263759/27., 15.46619297/27., SQRT3/8., -9.693949466/27.,
	0, 0, .25, -19.96676927/27., 22.64018827/27., -6.257417473/27.,
	0, 0, 0, .125, -13.39411601/27., 22.21162861/27.,
	0, 0, 0, 0, 1./16., -8.395206501/27.,
	0, 0, 0, 0, 0, 1./32.,
  };

  static const double dgp5_into_dgp5_refined_1[] =
  {
	1., SQRT3/2., 0, -8.929410675/27., 0, 5.596804334/27.,
	0, .5, 26.14263759/27., 15.46619297/27., -SQRT3/8., -9.693949466/27.,
	0, 0, .25, 19.96676927/27., 22.64018827/27., 6.257417473/27.,
	0, 0, 0, .125, 13.39411601/27., 22.21162861/27.,
	0, 0, 0, 0, 1./16., 8.395206501/27.,
	0, 0, 0, 0, 0, 1./32.,
  };

  static const double dgp6_into_dgp6_refined_0[] =
  {
	1., -SQRT3/2., 0, 8.929410675/27., 0, -5.596804334/27., 0.,
	0, .5, -26.14263759/27., 15.46619297/27., SQRT3/8., -9.693949466/27., -2.634608531/27.,
	0, 0, .25, -19.96676927/27., 22.64018827/27., -6.257417473/27., -10.20379496/27.,
	0, 0, 0, .125, -13.39411601/27., 22.21162861/27., -16.09772402/27.,
	0, 0, 0, 0, 1./16., -8.395206501/27., 18.25310333/27.,
	0, 0, 0, 0, 0, 1./32., -5.044891251/27.,
	0, 0, 0, 0, 0, 0, 1./64.,
  };

  static const double dgp6_into_dgp6_refined_1[] =
  {
	1., SQRT3/2., 0, -8.929410675/27., 0, 5.596804334/27., 0.,
	0, .5, 26.14263759/27., 15.46619297/27., -SQRT3/8., -9.693949466/27., 2.634608531/27.,
	0, 0, .25, 19.96676927/27., 22.64018827/27., 6.257417473/27., -10.20379496/27.,
	0, 0, 0, .125, 13.39411601/27., 22.21162861/27., 16.09772402/27.,
	0, 0, 0, 0, 1./16., 8.395206501/27., 18.25310333/27.,
	0, 0, 0, 0, 0, 1./32., 5.044891251/27.,
	0, 0, 0, 0, 0., 0., 1./64.,
  };
}


template <>
const double * const
FE_DGP<1>::Matrices::embedding[][GeometryInfo<1>::children_per_cell] =
{
  {FE_DGP_1d::dgp0_into_dgp0_refined_0, FE_DGP_1d::dgp0_into_dgp0_refined_1},
  {FE_DGP_1d::dgp1_into_dgp1_refined_0, FE_DGP_1d::dgp1_into_dgp1_refined_1},
  {FE_DGP_1d::dgp2_into_dgp2_refined_0, FE_DGP_1d::dgp2_into_dgp2_refined_1},
  {FE_DGP_1d::dgp3_into_dgp3_refined_0, FE_DGP_1d::dgp3_into_dgp3_refined_1},
  {FE_DGP_1d::dgp4_into_dgp4_refined_0, FE_DGP_1d::dgp4_into_dgp4_refined_1},
  {FE_DGP_1d::dgp5_into_dgp5_refined_0, FE_DGP_1d::dgp5_into_dgp5_refined_1},
  {FE_DGP_1d::dgp6_into_dgp6_refined_0, FE_DGP_1d::dgp6_into_dgp6_refined_1}
};

template <>
const unsigned int FE_DGP<1>::Matrices::n_embedding_matrices
= sizeof(FE_DGP<1>::Matrices::embedding) /
sizeof(FE_DGP<1>::Matrices::embedding[0]);


/*
 * These elements are not defined by interpolation, therefore there
 * are no interpolation matrices.
 */

template <>
const double * const FE_DGP<1>::Matrices::projection_matrices[][GeometryInfo<1>::children_per_cell] = {{0}};

template <>
const unsigned int FE_DGP<1>::Matrices::n_projection_matrices
= 0;


#else // #if deal_II_dimension
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // #if deal_II_dimension
