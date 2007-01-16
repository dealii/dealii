//----------------------------  injection_q_hierarchical.cc  ---------------------------
//    $Id: injection_q_hierarchical.cc 12732 2006-03-28 23:15:45Z wolf $
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  injection_q_hierarchical.cc  ---------------------------


char logname[] = "injection_q_hierarchical/output";


#include "../fe/injection_common.h"


template <int dim>
void test ()
{
  for (unsigned int i=1; i<4; ++i)
    for (unsigned int j=i; j<4; ++j)
      do_check (FE_Q_Hierarchical<dim>(i), FE_Q_Hierarchical<dim>(j));
}
