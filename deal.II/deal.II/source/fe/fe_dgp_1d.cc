//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
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
	1., 0.,
	0., .5,
  };

  static const double dgp1_into_dgp1_refined_1[] =
  {
	1., .5,
	0., .5,
  };

  static const double dgp2_into_dgp2_refined_0[] =
  {
	1., 0., -.375,
	0., .5, 0.,
	0., 0., .25,
  };

  static const double dgp2_into_dgp2_refined_1[] =
  {
	1., .5, 0.,
	0., .5, .75,
	0., 0., .25,
  };

  static const double dgp3_into_dgp3_refined_0[] =
  {
	1., 0., -.375, 0.,
	0., .5, 0., -15.1875/27.,
	0., 0., .25, 0.,
	0., 0., 0., .125,
  };

  static const double dgp3_into_dgp3_refined_1[] =
  {
	1., .5, 0., -.125,
	0., .5, .75, .375,
	0., 0., .25, 16.875/27.,
	0., 0., 0., .125,
  };

  static const double dgp4_into_dgp4_refined_0[] =
  {
	1., 0., -.375, 0., 3.1640625/27.,
	0., .5, 0., -15.1875/27., 0.,
	0., 0., .25, 0., -12.65625/27.,
	0., 0., 0., .125, 0.,
	0., 0., 0., 0., 1.6875/27.,
  };

  static const double dgp4_into_dgp4_refined_1[] =
  {
	1., .5, 0., -.125, 0.,
	0., .5, .75, .375, -.125,
	0., 0., .25, 16.875/27., 16.875/27.,
	0., 0., 0., .125, 11.8125/27.,
	0., 0., 0., 0., 1.6875/27.,
  };

  static const double dgp5_into_dgp5_refined_0[] =
  {
	1., 0., -.375, 0., 3.1640625/27., 0.,
	0., .5, 0., -15.1875/27., 0., 10.44140625/27.,
	0., 0., .25, 0., -12.65625/27., 0.,
	0., 0., 0., .125, 0., -8.859375/27.,
	0., 0., 0., 0., 1.6875/27., 0.,
	0., 0., 0., 0., 0., 0.84375/27.,
  };

  static const double dgp5_into_dgp5_refined_1[] =
  {
	1., .5, 0., -.125, 0., 1.6875/27.,
	0., .5, .75, .375, -.125, -5.0625/27.,
	0., 1.943780239e-09/27., 6.750000001/27., 16.875/27., 16.875/27., 4.21875/27.,
	0., -1.298171318e-09/27., 0., .125, 11.8125/27., 17.71875/27.,
	0., 0., 0., 0., 1.6875/27., 7.59375/27.,
	0., 0., 0., 0., 0., 0.84375/27.,
  };

  static const double dgp6_into_dgp6_refined_0[] =
  {
	1., 0., -10.12499999/27., 7.241569395e-09/27., 3.164062498/27., -5.02411607e-09/27., 0.553710937/27.,
	0., .5, -2.01244638e-08/27., -15.18750002/27., 5.331617674e-09/27., 10.44140626/27., 1.305959489e-09/27.,
	0., 0., 6.750000023/27., 2.013347609e-08/27., -12.65625001/27., -1.396873349e-08/27., 13.44726562/27.,
	0., 0., -1.701242498e-08/27., 3.374999985/27., 4.508076571e-09/27., -8.859374989/27., 1.102114197e-09/27.,
	0., 0., 8.95263339e-09/27., 7.969774316e-09/27., 1.687499998/27., -5.529893409e-09/27., -5.695312501/27.,
	0., 0., -3.09774831e-09/27., -2.756474137e-09/27., 0., 0.8437500019/27., 0.,
	0., 0., 0., 0., 0., 0., 0.421875/27.,
  };

  static const double dgp6_into_dgp6_refined_1[] =
  {
	1., 13.50000001/27., -6.32567268e-09/27., -3.37500001/27., 9.901835418e-09/27., 1.687500006/27., 6.235854535e-09/27.,
	0., 13.49999998/27., 20.25000002/27., 10.12500003/27., -3.375000025/27., -5.062500014/27., 1.265624985/27.,
	0., 2.037605333e-08/27., 6.749999982/27., 16.87499997/27., 16.87500003/27., 4.218750016/27., -6.328124983/27.,
	0., -1.534920617e-08/27., 1.32272523e-08/27., 3.375000021/27., 11.81249998/27., 17.71874999/27., 11.81249999/27.,
	0., 8.088959335e-09/27., -6.95629168e-09/27., -1.122813704e-08/27., 1.687500011/27., 7.593750006/27., 15.18750001/27.,
	0., -2.804903374e-09/27., 2.404689106e-09/27., 3.885094046e-09/27., -3.773966236e-09/27., 0.8437499978/27., 4.640624998/27.,
	0., 0., 0., 0., 0., 0., 0.4218750004/27.,
  };


};


template <>
const double * const
FE_DGP<1>::Matrices::embedding[][GeometryInfo<1>::children_per_cell] =
{
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
const double * const FE_DGP<1>::Matrices::projection_matrices[][GeometryInfo<1>::children_per_cell];

template <>
const unsigned int FE_DGP<1>::Matrices::n_projection_matrices
= 0;


#else // #if deal_II_dimension
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {}; };
#endif // #if deal_II_dimension
