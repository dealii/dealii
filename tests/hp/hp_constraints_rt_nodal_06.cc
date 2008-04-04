//----------------------------  hp_constraints_rt_nodal_06.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_constraints_rt_nodal_06.cc  ---------------------------


// check that computation of hp constraints works for RT elements correctly

char logname[] = "hp_constraints_rt_nodal_06/output";


#include "../hp/hp_constraints_common.h"
#include <fe/fe_raviart_thomas.h>


template <int dim>
void test ()
{
  if (dim == 1)
    return;
  
  hp::FECollection<dim> fe;
  std::vector<unsigned int> degrees;
  for (unsigned int i=0; i<7-dim; ++i)
    {
      fe.push_back (FE_RaviartThomasNodal<dim>(i));
      degrees.push_back (i);
    }
  
  test_interpolation  (fe, degrees);
}
