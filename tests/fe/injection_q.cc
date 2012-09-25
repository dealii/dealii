//----------------------------  injection_q.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2007, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  injection_q.cc  ---------------------------


char logname[] = "injection_q/output";


#include "injection_common.h"


template <int dim>
void test ()
{
  deallog << std::setprecision (4);
  for (unsigned int i=1; i<4; ++i)
    for (unsigned int j=i; j<4; ++j)
      do_check (FE_Q<dim>(i), FE_Q<dim>(j));
}
