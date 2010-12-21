//----------------------------  project_abf_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  project_abf_01.cc  ---------------------------


// check that VectorTools::project works for ABF elements correctly

char logname[] = "project_abf_01/output";


#include "../deal.II/project_common.h"


template <int dim>
void test ()
{
  if (dim > 1)
    for (unsigned int p=0; p<6-dim; ++p)
      test_no_hanging_nodes (FE_ABF<dim>(p), p);
}
