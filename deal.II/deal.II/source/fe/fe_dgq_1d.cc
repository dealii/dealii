//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------


// only compile in 1d
#if deal_II_dimension == 1

#include <fe/fe_dgq.h>


// Transfer matrices for finite elements

namespace FE_DGQ_1d
{
  static const double dgq0_into_dgq0_refined[] =
  {
	1., 1.
  };

  static const double dgq1_into_dgq1_refined[] =
  {
	1., 0.,
	13.5/27., 13.5/27.,
  };

  static const double dgq2_into_dgq2_refined[] =
  {
	1., 0., 0.,
	10.125/27., 20.25/27., -3.375/27.,
	0., 1., 0.,
  };

  static const double dgq3_into_dgq3_refined[] =
  {
	1., 0., 0., 0.,
	8.4375/27., 25.3125/27., -8.4375/27., 1.6875/27.,
	0., 1., 0., 0.,
	-1.6875/27., 15.1875/27., 15.1875/27., -1.6875/27.,
  };

  static const double dgq4_into_dgq4_refined[] =
  {
	1., 0., 0., 0., 0.,
	7.3828125/27., 29.53125/27., -14.765625/27., 5.90625/27., -1.0546875/27.,
	0., 1., 0., 0., 0.,
	-1.0546875/27., 12.65625/27., 18.984375/27., -4.21875/27., 0.6328125/27.,
	0., -2.91892343e-13/27., 1., 0., 0.,
  };


  static const double dgq0_refined_onto_dgq0[] =
  {
	0.5,
  };

  static const double dgq1_refined_onto_dgq1[] =
  {
	0.75, 0.5,
	-0.25, 0,
  };

  static const double dgq2_refined_onto_dgq2[] =
  {
	0.6875, 0.75, -0.1875,
	-0.09375, 0.375, 0.21875,
	0.1875, -0.25, -0.1875,
  };
  static const double dgq3_refined_onto_dgq3[] =
  {
	0.6875, 0.9375, -0.75, -0.0625,
	-0.055555556, 0.24305556, 0.61805556, 0.12268519,
	0.034722222, -0.11805556, 0.069444444, 0.085648148,
	-0.125, 0.1875, 0.1875, -0.0625,
  };

  static const double dgq4_refined_onto_dgq4[] =
  {
	0.72569444, 0.97222222, -1.0416667, 0.13888889, 0.017361111,
	-0.068223741, 0.30056424, 0.52115885, 0.28103299, 0.031873915,
	0.01953125, -0.052083333, -0.078125, 0.46875, 0.14192708,
	-0.0074598524, 0.036675347, -0.030924479, -0.066189236, 0.0014919705,
	0.086805556, -0.13888889, -0.20833333, 0.36111111, 0.086805556,
  };

  static const double dgq5_refined_onto_dgq5[] =
  {
	0.76529948, 0.96028646, -1.2858073, 0.42317708, 0.18717448, 0.10611979,
	-0.064940104, 0.29386458, 0.53121354, 0.25160417, 0.056143229, -0.012855417,
	0.036778646, -0.14359896, 0.12031771, 0.22294792, 0.53740365, 0.10436104,
	-0.0056692708, 0.022895833, 0.012755208, -0.16279167, 0.17783073, 0.076769167,
	-0.0047838542, 0.0020989583, 0.018526042, 0.016885417, -0.070075521, -0.017681042,
	-0.071940104, 0.16276042, 0.081380208, -0.32552083, -0.040690104, 0.037760417,
  };

  static const double dgq6_refined_onto_dgq6[] =
  {
	0.79257813, 0.984375, -1.6734375, 1.09375, -0.17226562, 0.196875, -0.065625,
	-0.056868389, 0.26796875, 0.56824544, 0.18953832, 0.14839681, -0.089322917, -0.019566816,
	0.03079829, -0.109375, 0.018229167, 0.41571502, 0.32044271, 0.328125, 0.038888889,
	-0.01965332, 0.0984375, -0.16918945, 0.109375, -0.1307373, 0.4921875, 0.11958008,
	-0.0025543338, -0.003125, 0.025520833, 0.0084619342, -0.081119792, -0.0072916667, 0.017283951,
	0.010613285, -0.02734375, 0.00079752604, 0.016814558, 0.058162435, -0.04921875, -0.018216508,
	0.062890625, -0.196875, 0.0984375, 0.13125, 0.12304687, -0.309375, -0.065625,
  };

  static const double dgq7_refined_onto_dgq7[] =
  {
	0.81508102, 1.0095262, -2.1220573, 2.1159252, -1.1456019, 0.72957682, -0.50176794, -0.037400897,
	-0.048702752, 0.24752223, 0.56041891, 0.24320734, 0.11116474, -0.063333627, -0.094508404, -0.012366406,
	0.028512974, -0.11272486, 0.058287533, 0.34450213, 0.37238149, 0.28423557, 0.10277054, 0.00045127394,
	-0.01556376, 0.056866734, 0.0023941134, -0.26908091, 0.36033153, 0.082017998, 0.54686868, 0.092819904,
	0.0063367825, -0.046119071, 0.12266792, -0.15184748, 0.13891781, -0.23206465, 0.23358346, 0.071870952,
	0.0077008445, -0.01866856, -0.01883456, 0.072907425, -0.023311084, -0.0010966409, -0.086150953, -0.010963118,
	-0.012439378, 0.040990099, -0.017439155, -0.050376157, 0.0042995099, 0.054999345, 0.038394855, -0.001831157,
	-0.050934606, 0.18647931, -0.15057292, -0.1043873, -0.016304977, 0.28934245, 0.015028935, -0.031932147,
  };
};



template <>
const double * const FE_DGQ<1>::Matrices::embedding[] =
{	
      FE_DGQ_1d::dgq0_into_dgq0_refined,
      FE_DGQ_1d::dgq1_into_dgq1_refined,
      FE_DGQ_1d::dgq2_into_dgq2_refined,
      FE_DGQ_1d::dgq3_into_dgq3_refined,
      FE_DGQ_1d::dgq4_into_dgq4_refined
};



template <>
const unsigned int FE_DGQ<1>::Matrices::n_embedding_matrices
  = sizeof(FE_DGQ<1>::Matrices::embedding) /
    sizeof(FE_DGQ<1>::Matrices::embedding[0]);




template <>
const double * const FE_DGQ<1>::Matrices::projection_matrices[] =
{
      FE_DGQ_1d::dgq0_refined_onto_dgq0,
      FE_DGQ_1d::dgq1_refined_onto_dgq1,
      FE_DGQ_1d::dgq2_refined_onto_dgq2,
      FE_DGQ_1d::dgq3_refined_onto_dgq3,
      FE_DGQ_1d::dgq4_refined_onto_dgq4,
      FE_DGQ_1d::dgq5_refined_onto_dgq5,
      FE_DGQ_1d::dgq6_refined_onto_dgq6,
      FE_DGQ_1d::dgq7_refined_onto_dgq7
};



template <>
const unsigned int FE_DGQ<1>::Matrices::n_projection_matrices
  = sizeof(FE_DGQ<1>::Matrices::projection_matrices) /
    sizeof(FE_DGQ<1>::Matrices::projection_matrices[0]);


#else // #if deal_II_dimension
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {}; };
#endif // #if deal_II_dimension == 1
