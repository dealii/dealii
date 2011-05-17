//----------------------------  hp_constraints_rt_04.cc  ---------------------------
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
//----------------------------  hp_constraints_rt_04.cc  ---------------------------


// check that computation of hp constraints works for RT elements correctly

char logname[] = "hp_constraints_rt_04/output";


#include "../hp/hp_constraints_common.h"
#include <deal.II/fe/fe_raviart_thomas.h>


template <int dim>
void test ()
{
  if (dim == 1)
    return;
  
  hp::FECollection<dim> fe;
  for (unsigned int i=1; i<4; ++i)
    fe.push_back (FE_RaviartThomas<dim>(i));
  test_with_2d_deformed_mesh  (fe);
}
