//----------------------------  injection_dgp_nonparametric.cc  ---------------------------
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
//----------------------------  injection_dgp_nonparametric.cc  ---------------------------


char logname[] = "injection_dgp_nonparametric/output";


#include "../fe/injection_common.h"


template <int dim>
void test ()
{
  for (unsigned int i=1; i<4; ++i)
    for (unsigned int j=i; j<4; ++j)
      do_check (FE_DGPNonparametric<dim>(i), FE_DGPNonparametric<dim>(j));
}
