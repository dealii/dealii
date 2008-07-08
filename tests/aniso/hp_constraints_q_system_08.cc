//----------------------------  hp_constraints_q_system_08.cc  ---------------------------
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
//----------------------------  hp_constraints_q_system_08.cc  ---------------------------


// check that computation of hp constraints works for FESystem(FE_Q) elements correctly
// on a uniformly refined mesh for functions of degree q

char logname[] = "hp_constraints_q_system_08/output";


#include "../hp/hp_constraints_common.h"


template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  for (unsigned int i=1; i<4; ++i)
    fe.push_back (FESystem<dim>(FE_Q<dim>(i), 1,
				FE_DGQ<dim>(i+1), 1));
  
  test_with_hanging_nodes_random_aniso (fe);
}
