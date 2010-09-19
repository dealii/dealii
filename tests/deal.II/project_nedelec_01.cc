//----------------------------  project_nedelec_01.cc  ---------------------------
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
//----------------------------  project_nedelec_01.cc  ---------------------------


// check that VectorTools::project works for Nedelec elements correctly

char logname[] = "project_nedelec_01/output";


#include "project_common.h"


template <int dim>
void test ()
{
  if (dim > 1)
				     // only p=1 implemented at present
    for (unsigned int p=1; p<2; ++p)
      test_no_hanging_nodes (FE_Nedelec<dim>(p-1), p, 1);
}
