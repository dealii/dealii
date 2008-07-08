//----------------------------  hp_constraints_dgp_monomial_08.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_constraints_dgp_monomial_08.cc  ---------------------------


// check that computation of hp constraints works for DGP_Monomial elements correctly

char logname[] = "hp_constraints_dgp_monomial_08/output";


#include "../hp/hp_constraints_common.h"


template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  for (unsigned int i=0; i<4; ++i)
    fe.push_back (FE_DGPMonomial<dim>(i));
  
  test_with_hanging_nodes_random_aniso (fe);
}
