//----------------------------  hp_constraints_q_system_x_01.cc  ---------------------------
//    $Id: hp_constraints_q_system_x_01.cc 12732 2006-03-28 23:15:45Z wolf $
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_constraints_q_system_x_01.cc  ---------------------------


// check that computation of hp constraints works for FESystem(FE_Q) elements correctly
// on a uniformly refined mesh for functions of degree q

// these tests check that we can deal with FESystem(FE_Q(p),FE_DGQ(q)) for
// different p,q

char logname[] = "hp_constraints_q_system_x_01/output";


#include "hp_constraints_common.h"


template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  for (unsigned int i=1; i<4; ++i)
    for (unsigned int j=0; j<4; ++j)
      fe.push_back (FESystem<dim>(FE_Q<dim>(i), 1,
				  FE_DGQ<dim>(j), 1));
  
  test_no_hanging_nodes (fe);
}
