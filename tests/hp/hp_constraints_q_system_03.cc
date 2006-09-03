//----------------------------  hp_constraints_q_system_03.cc  ---------------------------
//    $Id: hp_constraints_q_system_03.cc 12732 2006-03-28 23:15:45Z wolf $
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_constraints_q_system_03.cc  ---------------------------


// check that computation of hp constraints works for FESystem(FE_Q) elements correctly
// on a uniformly refined mesh for functions of degree q

char logname[] = "hp_constraints_q_system_03/output";


#include "hp_constraints_common.h"


template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  for (unsigned int i=1; i<4; ++i)
    fe.push_back (FESystem<dim>(FE_Q<dim>(i), 1,
				FE_DGQ<dim>(i+1), 1));
  
  test_with_wrong_face_orientation (fe);
}
