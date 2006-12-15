//----------------------------  hp_constraints_q_system_06.cc  ---------------------------
//    $Id: hp_constraints_q_system_06.cc 12732 2006-03-28 23:15:45Z wolf $
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_constraints_q_system_06.cc  ---------------------------


// check that computation of hp constraints works for FESystem(FE_Q) elements correctly
// on a uniformly refined mesh for functions of degree q

char logname[] = "hp_constraints_q_system_06/output";


#include "hp_constraints_common.h"


template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  std::vector<unsigned int> degrees;
  for (unsigned int i=1; i<7-dim; ++i)
    {
				       // in contrast to all the other
				       // hp_constraints_q_system_* tests, we
				       // use two continuous elements here,
				       // since we can't interpolate DGQ
				       // elements
      fe.push_back (FESystem<dim>(FE_Q<dim>(i), 1,
				  FE_Q<dim>(i+1), 1));
      degrees.push_back (i);
    }
  
  test_interpolation  (fe, degrees);
}
