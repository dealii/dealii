//----------------------------  hp_constraints_q_06.cc  ---------------------------
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
//----------------------------  hp_constraints_q_06.cc  ---------------------------


// check that computation of hp constraints works for Q elements correctly

char logname[] = "hp_constraints_q_06/output";


#include "hp_constraints_common.h"


template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  std::vector<unsigned int> degrees;
  for (unsigned int i=1; i<7-dim; ++i)
    {
      fe.push_back (FE_Q<dim>(i));
      degrees.push_back (i);
    }
  
  test_interpolation  (fe, degrees);
}
