//----------------------------  hp_constraints_q_system_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_constraints_q_system_01.cc  ---------------------------


// this is reduced from hp_constraints_q_system_x_01, where we test that we
// can deal with FESystem(FE_Q(p),FE_DGQ(q)) for different p,q. note that for
// fixed p but varying q, neither of the two elements will dominate the other

char logname[] = "crash_13/output";


#include "hp_constraints_common.h"


template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  fe.push_back (FESystem<dim>(FE_Q<dim>(1), 1,
			      FE_DGQ<dim>(0), 1));
  fe.push_back (FESystem<dim>(FE_Q<dim>(1), 1,
			      FE_DGQ<dim>(1), 1));
  
  test_no_hanging_nodes (fe);
}
