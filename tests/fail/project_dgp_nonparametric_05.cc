//----------------------------  project_dgp_nonparametric_05.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  project_dgp_nonparametric_05.cc  ---------------------------


// check that VectorTools::project works for DGPNonparametric elements
// correctly on a uniformly refined mesh for functions of degree q

char logname[] = "project_dgp_nonparametric_05/output";


#include "../deal.II/project_common.h"


template <int dim>
void test ()
{
  for (unsigned int p=0; p<6-dim; ++p)
    test_with_2d_deformed_refined_mesh (FE_DGPNonparametric<dim>(p), p);
}
