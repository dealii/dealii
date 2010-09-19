//----------------------------  project_dgq_04.cc  ---------------------------
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
//----------------------------  project_dgq_04.cc  ---------------------------


// check that VectorTools::project works for DGQ elements correctly

char logname[] = "project_dgq_04/output";


#include "project_common.h"


template <int dim>
void test ()
{
  for (unsigned int p=0; p<6-dim; ++p)
    test_with_2d_deformed_mesh (FE_DGQ<dim>(p), p);
}
