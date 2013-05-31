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


char logname[] = "injection_q_iso_q1/output";


#include "injection_common.h"
#include <deal.II/fe/fe_q_iso_q1.h>

template <int dim>
void test ()
{
  deallog << std::setprecision (4);
  for (unsigned int i=1; i<4; ++i)
    for (unsigned int j=i; j<4; ++j)
      if (j%i == 0)
        do_check (FE_Q_iso_Q1<dim>(i), FE_Q_iso_Q1<dim>(j));
}
