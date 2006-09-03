//----------------------------  hp_constraints_q_04.cc  ---------------------------
//    $Id: hp_constraints_q_04.cc 12732 2006-03-28 23:15:45Z wolf $
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_constraints_q_04.cc  ---------------------------


// check that computation of hp constraints works for Q elements correctly

char logname[] = "hp_constraints_q_04/output";


#include "hp_constraints_common.h"


template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  for (unsigned int i=1; i<4; ++i)
    fe.push_back (FE_Q<dim>(i));
  test_with_2d_deformed_mesh  (fe);
}
