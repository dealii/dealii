//----------------------------  hp_constraints_dgp_nonparametric_06.cc  ---------------------------
//    $Id: hp_constraints_dgp_nonparametric_06.cc 12732 2006-03-28 23:15:45Z wolf $
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_constraints_dgp_nonparametric_06.cc  ---------------------------


// check that computation of hp constraints works for DGPNonparametric elements
// correctly on a uniformly refined mesh for functions of degree q

char logname[] = "hp_constraints_dgp_nonparametric_06/output";


#include "hp_constraints_common.h"


template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  std::vector<unsigned int> degrees;
  for (unsigned int i=0; i<7-dim; ++i)
    {
      fe.push_back (FE_DGPNonparametric<dim>(i));
      degrees.push_back (i);
    }
  
  test_interpolation  (fe, degrees);
}
